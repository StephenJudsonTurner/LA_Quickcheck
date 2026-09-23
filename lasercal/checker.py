"""Automatic first-pass reduction ("data checker").

Given either a folder of Thermo Qtegra per-analysis CSVs or a MATLAB-style
data folder (data workbook + RunOrder.xlsx), choose reduction parameters
from the data itself, run the reduction, and write the standard plots and an
HTML report.  Every choice is written to the output folder in the same XLSX
files the tool normally reads, so a second pass can edit them by hand.

Rules used for the automatic choices (see choose_* functions):

* windows       laser-on / laser-off from the 25Mg trace of each analysis;
                background = [1 s, laser-on - 2 s]; signal starts 2 s after
                laser-on (the first second carries a surface Pb/Cu spike and
                the sweep-timing artefact) and lasts the standards' median
                ablation length - 4 s, the same for every sample. An unknown
                whose ablation is shorter gets its own (ablation - 2 s) window
                and is flagged (down-hole mismatch). A second background only
                if >= 12 s of clean gas blank follow the ablation.
* drift         pooled: every recognised standard in every bracket, each
                analysis weighted by its counting precision, the element's
                ratio to its normaliser divided by that standard's run mean;
                pchip through the bracket means when they are good to 1.5 %,
                otherwise a quadratic (core.DRIFT_SE_SWITCH_PCT).
* calibration   normalise to Al (Ca for 27Al) when Al was measured, else Ca;
                slope = inverse-variance weighted mean of log(target/measured)
                over the BHV/BCR/BIR replicates (each standard counts by its
                counting precision, not by its concentration); replicates
                with > 10 % counting error are dropped; GeoRem values.  A
                calibrant replicate that departs from its own standard's other
                replicates by > 4 robust sigma and > 10 % is excluded too
                (never more than a third of a standard's points) and listed.
* reporting     the counting-statistics floor is printed next to every RSD,
                bias is shown only where that floor is < 5 %, and the
                apparent sensitivity of every standard relative to the
                calibration (reference-value / matrix consistency) and the
                down-hole fractionation slopes of the standards are tabulated.
"""
from __future__ import annotations

import html
import os
import re
import shutil
import sys
from datetime import datetime
from typing import Callable, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

from . import core, io_xlsx, qtegra
from .names import clean_sample_names
from .session import ReductionSession, UIState, IntervalsCreated

Log = Callable[[str], None]

SETTINGS_FILES = ('RunOrder.xlsx', 'Intervals.xlsx', 'DriftSelections.xlsx', 'CalibrationSelections.xlsx',
                  'SampleWindows.xlsx', 'DriftCorrected_All.xlsx', 'ReducedDataExport.xlsx', 'SourceFiles.xlsx',
                  'LAICPMS_Standard_Values.xlsx', 'DwellTimes.xlsx', 'StandardConsistency.xlsx', 'DownholeSlopes.xlsx',
                  'NormaliserChoice.xlsx')

SIGNAL_DELAY_S = 2.0     # signal starts this long after laser-on
SIGNAL_TAIL_S = 2.0      # and ends this long before laser-off (of the standards' median ablation)
BIAS_FLOOR_PCT = 5.0     # bias is reported only where the counting floor is below this
NORM_SWITCH_MARGIN_PCT = 0.2   # a normaliser must beat the default's score by this to be chosen
CONSISTENCY_FLAG_PCT = 5.0


def _log_default(msg: str):
    print(msg, flush=True)


# ---------------------------------------------------------------- input detection
def detect_input(path: str) -> Tuple[str, str]:
    """-> ('qtegra', folder) | ('xlsx', data workbook path)."""
    if os.path.isdir(path):
        csvs = qtegra.list_files(path)
        if csvs:
            return 'qtegra', path
        xl = [f for f in os.listdir(path) if f.lower().endswith('.xlsx') and f not in SETTINGS_FILES]
        if len(xl) == 1 and os.path.isfile(os.path.join(path, 'RunOrder.xlsx')):
            return 'xlsx', os.path.join(path, xl[0])
        raise ValueError('Folder must contain Qtegra CSV files, or exactly one data workbook plus RunOrder.xlsx.')
    if path.lower().endswith('.xlsx'):
        return 'xlsx', path
    raise ValueError(f'Unrecognised input: {path}')


# ---------------------------------------------------------------- choices
def choose_windows(windows: List[pd.DataFrame], thresh: float = core.START_THRESH,
                   iso: str = 'x25Mg', std_idx: Optional[np.ndarray] = None) -> Tuple[pd.DataFrame, Dict]:
    """Per-sample Intervals from the laser-on/off times in each window."""
    on_t, off_t, end_t = [], [], []
    for T in windows:
        t = T['Time'].to_numpy(float)
        mg = T[iso].to_numpy(float) if iso in T.columns else np.zeros_like(t)
        hit = np.flatnonzero(mg >= thresh)
        on_t.append(t[hit[0]] if hit.size else np.nan)
        off_t.append(t[hit[-1]] if hit.size else np.nan)
        end_t.append(t[-1])
    on_t, off_t, end_t = map(np.array, (on_t, off_t, end_t))
    ablation = off_t - on_t
    # The whole ablation minus its first SIGNAL_DELAY_S and last SIGNAL_TAIL_S seconds, sized on
    # the standards (the primaries are counting-limited: a 20 s cap threw away 30 % of the counts).
    ref = ablation[std_idx] if std_idx is not None and len(std_idx) else ablation
    sig_len = float(np.nanmedian(ref)) - SIGNAL_DELAY_S - SIGNAL_TAIL_S
    sig_len = max(sig_len, 4.0)
    own_len = np.maximum(ablation - SIGNAL_DELAY_S - SIGNAL_TAIL_S, 2.0)
    short = np.isfinite(ablation) & (own_len < sig_len - 1.0)
    length = np.where(short, own_len, sig_len)
    tail = end_t - off_t
    use_back2 = bool(np.nanmin(tail) >= 12.0)
    iv = pd.DataFrame(np.nan, index=range(len(windows)), columns=core.INTERVAL_COLUMNS)
    iv['back_start1'] = 1.0
    iv['back_stop1'] = on_t - 2.0
    iv['signal_start'] = on_t + SIGNAL_DELAY_S
    iv['signal_stop'] = on_t + SIGNAL_DELAY_S + length
    if use_back2:
        iv['back_start2'] = off_t + 4.0
        iv['back_stop2'] = end_t - 0.5
    info = {'laser_on_median_s': float(np.nanmedian(on_t)), 'ablation_median_s': float(np.nanmedian(ref)),
            'signal_length_s': sig_len, 'second_background': use_back2,
            'n_no_laser_on': int(np.isnan(on_t).sum()),
            'short_signal': [int(i) for i in np.flatnonzero(short)],
            'short_signal_len': [float(own_len[i]) for i in np.flatnonzero(short)]}
    return iv, info


def choose_drift(sets: Dict[str, np.ndarray]) -> Tuple[str, str, str]:
    """-> (standard, method, reason).  Pooled over every recognised standard."""
    present = [k for k in core.STANDARD_TAGS if sets.get(k, np.array([])).size > 0]
    if not present:
        return 'ALL', 'none', 'no standard found; no drift correction'
    n_groups = len(core.standard_brackets(sets))
    if n_groups < 2:
        return 'ALL', 'none', 'standards run in a single bracket; no drift correction possible'
    n_per = ', '.join(f'{k} x{sets[k].size}' for k in present)
    return 'ALL', 'auto', (f'pooled over {n_per} in {n_groups} brackets, each analysis weighted by its counting '
                           f'precision; pchip through the bracket means where they are good to '
                           f'{core.DRIFT_SE_SWITCH_PCT:g} %, quadratic otherwise')


def choose_calibration(sig_vars: List[str], norm_map) -> pd.DataFrame:
    have_al = norm_map['Al']['iso'] in sig_vars
    have_ca = norm_map['Ca']['iso'] in sig_vars
    rows = []
    for v in sig_vars:
        if have_al and v != norm_map['Al']['iso']:
            key = 'Al'
        elif have_ca and v != norm_map['Ca']['iso']:
            key = 'Ca'
        else:
            key = 'Al' if have_al else ('Ca' if have_ca else 'Si')
        rows.append({'NormElement': key, 'ForceInterceptZero': True, 'StandardSet': 'GeoRem', 'CalExclusions': 'None',
                     'Weighting': 'ivw', 'Calibrants': ' '.join(core.default_calibrants(v))})
    return pd.DataFrame(rows, index=sig_vars, dtype=object)


def normaliser_scores(s: ReductionSession) -> pd.DataFrame:
    """For the normaliser the session is currently set to: per analyte, the RMS
    over standards of the replicate RSD of the drift-corrected ratio (precision)
    and the spread (std of log) of the apparent sensitivity across the standards
    with reference values (accuracy), both in %, using only standards whose
    counting error for the ratio is below BIAS_FLOOR_PCT."""
    ratio = s.measured_ratios(); err = s.ratio_error_pct(); std = s.std_vals
    rows = {}
    for j, el in enumerate(s.sig_vars):
        norm_iso = s.sig_vars[s.norm_idx[j]]
        rsds, sens = [], []
        for k in core.STANDARD_TAGS:
            idx = s.sets.get(k, np.array([], int))
            if idx.size < 3 or np.nanmedian(err[idx, j]) > BIAS_FLOOR_PCT:
                continue
            v = ratio[idx, j]
            if not np.all(np.isfinite(v)) or np.mean(v) <= 0:
                continue
            rsds.append(100 * np.std(v, ddof=1) / np.mean(v))
            ri = core.std_row_index(std, k, 'GeoRem')
            if ri is None or el not in std.columns or norm_iso not in std.columns:
                continue
            tgt = pd.to_numeric(pd.Series([std[el].iloc[ri], std[norm_iso].iloc[ri]]), errors='coerce').to_numpy(float)
            if np.all(np.isfinite(tgt)) and tgt[0] > 0 and tgt[1] > 0:
                sens.append(np.log(np.mean(v) / (tgt[0] / tgt[1])))
        prec = float(np.sqrt(np.mean(np.square(rsds)))) if rsds else np.nan
        acc = float(100 * np.std(sens)) if len(sens) >= 3 else np.nan
        rows[el] = {'precision': prec, 'accuracy': acc,
                    'score': float(np.sqrt(np.nansum([prec ** 2, acc ** 2]))) if rsds else np.nan,
                    'n_std_prec': len(rsds), 'n_std_acc': len(sens)}
    return pd.DataFrame(rows).T


def choose_normalisers(s: ReductionSession, cal: pd.DataFrame, log: Log) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Try every available normaliser on the whole run and keep, per analyte, the
    one with the lowest combined score; the default (from choose_calibration)
    is kept unless another beats it by NORM_SWITCH_MARGIN_PCT.  Returns the
    updated calibration preferences and the score table."""
    cands = list(s.show_choices)
    default = cal['NormElement'].copy()
    scores = {}
    for key in cands:
        trial = cal.copy()
        for el in trial.index:
            trial.loc[el, 'NormElement'] = key if s.norm_map[key]['iso'] != el else next(
                (k for k in cands if s.norm_map[k]['iso'] != el), key)
        io_xlsx.write_calibration_selections(s.cal_prefs_path, trial)
        s.run_drift(); s.calibrate()
        scores[key] = normaliser_scores(s)
    table = pd.concat({k: v[['precision', 'accuracy', 'score']] for k, v in scores.items()}, axis=1)
    chosen = {}
    for el in cal.index:
        d = default[el]
        best, best_score = d, scores[d].loc[el, 'score'] if d in scores else np.inf
        for key in cands:
            if s.norm_map[key]['iso'] == el:
                continue
            sc = scores[key].loc[el, 'score']
            if np.isfinite(sc) and (not np.isfinite(best_score) or sc < best_score - NORM_SWITCH_MARGIN_PCT):
                best, best_score = key, sc
        chosen[el] = best
    cal = cal.copy(); cal['NormElement'] = pd.Series(chosen)
    io_xlsx.write_calibration_selections(s.cal_prefs_path, cal)
    flat = table.copy(); flat.columns = [f'{a}_{b}' for a, b in table.columns]
    flat['chosen'] = pd.Series(chosen)
    io_xlsx.write_row_table(os.path.join(s.folder, 'NormaliserChoice.xlsx'), flat, 'Sheet1')
    changed = {el: (default[el], chosen[el]) for el in cal.index if chosen[el] != default[el]}
    log('normaliser per analyte (score = replicate RSD and cross-standard spread in quadrature, %): ' +
        ', '.join(f'{k}: {int((cal.NormElement == k).sum())}' for k in cands) +
        (('; changed from the default: ' + ', '.join(f'{el} {a}->{b}' for el, (a, b) in changed.items())) if changed else '; none changed'))
    return cal, flat


def auto_exclusions(s: ReductionSession, n_sigma: float = 4.0, min_rel: float = 0.10) -> Dict[str, List[str]]:
    """Gross outliers among the replicates of one calibrant standard, per analyte.

    A replicate is excluded when its measured/target ratio is more than
    ``n_sigma`` robust sigma (MAD) from the median of its own standard's
    replicates AND at least ``min_rel`` (10 %) away from it; never more than
    a third of a standard's points.  Systematic standard-to-standard offsets
    are deliberately NOT treated as outliers (they are reference-value issues
    for the user to judge).
    """
    out = {}
    for el in s.cal_elements:
        x, y, names = el.x_plot, el.y_plot, np.array(el.names_plot)
        chosen = []
        for g in ('bhv', 'bcr', 'bir'):
            m = np.array([g in nm.lower() for nm in names]) & np.isfinite(x) & np.isfinite(y) & (x > 0) & (y > 0)
            if m.sum() < 3:
                continue
            r = np.log(x[m] / y[m])                     # log(measured/target) per replicate
            med = np.median(r)
            mad = 1.4826 * np.median(np.abs(r - med))
            dev = np.abs(r - med)
            flagged = (dev > n_sigma * max(mad, 1e-12)) & (dev > np.log(1 + min_rel))
            k_max = int(m.sum() // 3)
            if flagged.sum() > k_max:
                keep = np.argsort(-dev)[:k_max]
                f2 = np.zeros_like(flagged); f2[keep] = True; flagged = f2 & flagged
            chosen += list(names[m][flagged])
        if chosen:
            out[el.analyte] = chosen
    return out


# ---------------------------------------------------------------- run
def run_checker(input_path: str, out_dir: Optional[str] = None, std_file: Optional[str] = None,
                log: Log = _log_default) -> str:
    kind, src = detect_input(input_path)
    if out_dir is None:
        base = src if kind == 'qtegra' else os.path.dirname(src)
        out_dir = os.path.join(os.path.dirname(base.rstrip('\\/')),
                               os.path.basename(base.rstrip('\\/')) + '_checked')
    data_dir = os.path.join(out_dir, 'data')
    plots_dir = os.path.join(out_dir, 'plots')
    os.makedirs(data_dir, exist_ok=True)
    for f in os.listdir(data_dir):
        if f.endswith('.xlsx'):
            os.remove(os.path.join(data_dir, f))
    if std_file is None:
        std_file = bundled_standards_path()
    log(f'input: {src}  ({kind})')
    log(f'output: {out_dir}')

    # ---- build the tool-format data folder
    if kind == 'qtegra':
        name = os.path.basename(src.rstrip('\\/'))
        info = qtegra.build_data_folder(src, data_dir, name)
        data_path = info['data']
        log(f'{info["n_samples"]} analyses read from CSV files')
    else:
        data_path = os.path.join(data_dir, os.path.basename(src))
        shutil.copy2(src, data_path)
        shutil.copy2(os.path.join(os.path.dirname(src), 'RunOrder.xlsx'), os.path.join(data_dir, 'RunOrder.xlsx'))
        sw = os.path.join(os.path.dirname(src), 'SampleWindows.xlsx')
        if os.path.isfile(sw):
            shutil.copy2(sw, data_dir)
        log('data workbook + RunOrder.xlsx copied')
    shutil.copy2(std_file, os.path.join(data_dir, 'LAICPMS_Standard_Values.xlsx'))

    # ---- first load just to see the windows (Intervals.xlsx gets created + we stop)
    s = ReductionSession(data_path, std_file=std_file, ui=UIState())
    try:
        s.load()
    except IntervalsCreated:
        pass
    # windows/sets are populated even when Intervals.xlsx was just created
    sets = s.sets
    for k in ('BHV', 'BCR', 'BIR'):
        if sets[k].size == 0:
            raise ValueError(f'No {k} analyses found (names must contain "{k.lower()}"); cannot calibrate.')
    log(f'{s.n_samples} samples, {len(s.sig_vars)} analytes; standards: ' +
        ', '.join(f'{k} x{sets[k].size}' for k in core.STANDARD_TAGS if sets[k].size))

    iv, winfo = choose_windows(s.windows, std_idx=core.all_standard_idx(sets))
    iv.index = s.display_names
    io_xlsx.write_intervals(s.intervals_path, iv)
    log(f'windows: laser-on at {winfo["laser_on_median_s"]:.1f} s, standards ablate {winfo["ablation_median_s"]:.1f} s; '
        f'signal = laser-on + {SIGNAL_DELAY_S:g} s for {winfo["signal_length_s"]:.1f} s; second background: {winfo["second_background"]}')
    if winfo['short_signal']:
        log(f'{len(winfo["short_signal"])} analyses ablate for less than that and get their own shorter window (flagged in the report)')

    dstd, dmeth, dreason = choose_drift(sets)
    io_xlsx.write_drift_selections(s.drift_prefs_path, pd.DataFrame(
        {'Standard': dstd, 'Interpolation': dmeth, 'Exclusions': ''}, index=s.sig_vars, dtype=object))
    log(f'drift: {dreason}')
    cal = choose_calibration(s.sig_vars, s.norm_map)
    io_xlsx.write_calibration_selections(s.cal_prefs_path, cal)
    special = '; '.join(f'{v}: {c}' for v, c in cal.Calibrants.items() if c != ' '.join(core.PRIMARY_STDS))
    log(f'calibration: normalised to {cal.NormElement.mode()[0]} (Ca for the Al isotope), inverse-variance weighted '
        f'log-mean sensitivity over BHV/BCR/BIR, GeoRem values, replicates with > {core.MAX_CAL_ERR_PCT:g} % counting error dropped'
        + (f'; calibrants overridden for {special}' if special else ''))

    # ---- pick the normaliser per analyte from the standards, then reduce, find outliers, reduce again
    s = ReductionSession(data_path, std_file=std_file, ui=UIState()).load()
    cal, _ = choose_normalisers(s, cal, log)
    s = ReductionSession(data_path, std_file=std_file, ui=UIState()).load()
    s.run_drift(); s.calibrate()
    excl = auto_exclusions(s)
    if excl:
        cal = s.load_or_init_cal_prefs()
        for el, names in excl.items():
            cal.loc[el, 'CalExclusions'] = ' '.join(names)
        io_xlsx.write_calibration_selections(s.cal_prefs_path, cal)
        log('auto-excluded calibrant points: ' + '; '.join(f'{el}: {" ".join(n)}' for el, n in excl.items()))
    else:
        log('no calibrant outliers beyond 4 robust sigma')
    return reduce_and_report(data_path, out_dir, std_file, log, kind, src, winfo, (dstd, dmeth, dreason), excl)


def rerun_folder(out_dir: str, std_file: Optional[str] = None, log: Log = _log_default) -> str:
    """Re-reduce a *_checked folder after its data/*.xlsx settings were edited by hand."""
    data_dir = os.path.join(out_dir, 'data')
    if not os.path.isdir(data_dir):
        raise ValueError(f'{out_dir} has no data/ sub-folder (pick a *_checked folder).')
    std_file = std_file or (os.path.join(data_dir, 'LAICPMS_Standard_Values.xlsx')
                            if os.path.isfile(os.path.join(data_dir, 'LAICPMS_Standard_Values.xlsx')) else bundled_standards_path())
    skip = set(SETTINGS_FILES)
    data = [f for f in os.listdir(data_dir) if f.endswith('.xlsx') and f not in skip]
    if len(data) != 1:
        raise ValueError(f'expected one data workbook in {data_dir}, found {data}')
    data_path = os.path.join(data_dir, data[0])
    log(f're-running {data_path} with the settings in data/*.xlsx')
    s = ReductionSession(data_path, std_file=std_file, ui=UIState()).load()
    dp = s.load_or_init_drift_prefs()
    drift = (dp.Standard.mode()[0], dp.Interpolation.mode()[0], 'as set in data/DriftSelections.xlsx')
    cal = s.load_or_init_cal_prefs()
    excl = {a: core.parse_cal_exclusions(v) for a, v in cal.CalExclusions.items() if core.parse_cal_exclusions(v)}
    iv = s.intervals
    lens = (iv['signal_stop'] - iv['signal_start']).to_numpy(float)
    med = float(np.nanmedian(lens))
    winfo = {'laser_on_median_s': float(np.nan), 'ablation_median_s': float(np.nan),
             'signal_length_s': med, 'second_background': bool(iv['back_start2'].notna().any()),
             'short_signal': [int(i) for i in np.flatnonzero(np.isfinite(lens) & (lens < med - 0.5))],
             'short_signal_len': [float(l) for l in lens[np.isfinite(lens) & (lens < med - 0.5)]]}
    return reduce_and_report(data_path, out_dir, std_file, log, 'rerun', data_path, winfo, drift, excl)


def reduce_and_report(data_path, out_dir, std_file, log, kind, src, winfo, drift, excl) -> str:
    plots_dir = os.path.join(out_dir, 'plots')
    s = ReductionSession(data_path, std_file=std_file, ui=UIState()).load()
    s.run_drift(); s.calibrate()
    s.export_drift_corrected()
    out_xlsx = s.export_calibration()
    log(f'wrote {out_xlsx}')

    # ---- plots
    import plot_calibrations as pc
    import plot_standard_recovery as pr
    cal_dir = os.path.join(plots_dir, 'calibration_curves'); os.makedirs(cal_dir, exist_ok=True)
    rec_dir = os.path.join(plots_dir, 'standard_recovery'); os.makedirs(rec_dir, exist_ok=True)
    for d in (cal_dir, rec_dir):
        for f in os.listdir(d):
            if f.endswith('.png') or f.endswith('.xlsx'):
                os.remove(os.path.join(d, f))
    err = core.compute_sem_backcorr_normalized(s.windows, s.intervals, s.sig_vars, s.norm_map,
                                               s.load_or_init_cal_prefs(), stat='sem')
    cal_rows, rec_rows, cal_pngs, rec_pngs = [], [], [], []
    for j, el in enumerate(s.cal_elements):
        p1 = os.path.join(cal_dir, f'{el.analyte}.png'); pc.plot_one(el, s, p1); cal_pngs.append(p1)
        p2 = os.path.join(rec_dir, f'{el.analyte}.png'); rec_rows += pr.plot_one(s, el, j, err, 2.0, 'SEM', p2); rec_pngs.append(p2)
        n_exc = len([n for n in el.names_plot if n in el.excluded])
        cal_rows.append({'Analyte': el.analyte, 'Normaliser': s.norm_map[el.norm_key]['iso'], 'StandardSet': el.set_name,
                         'Slope': el.slope, 'Intercept': el.intercept, 'R2': el.r2,
                         'nIncluded': len(el.names_plot) - n_exc, 'nExcluded': n_exc,
                         'Excluded': ' '.join(el.excluded) if el.excluded else 'None'})
    io_xlsx.write_table(os.path.join(cal_dir, 'calibration_stats.xlsx'), pd.DataFrame(cal_rows), 'Sheet1')
    io_xlsx.write_table(os.path.join(rec_dir, 'recovery_stats.xlsx'), pd.DataFrame(rec_rows), 'Sheet1')
    pc.contact_sheet(cal_pngs, os.path.join(cal_dir, 'index.png'))
    pc.contact_sheet(rec_pngs, os.path.join(rec_dir, 'index.png'))
    log(f'{len(cal_pngs)} calibration and {len(rec_pngs)} recovery plots written')

    # ---- counting floor, consistency of the standards, down-hole slopes
    data_dir = os.path.dirname(data_path)
    floor = pd.DataFrame(s.ratio_error_pct(), index=s.display_names, columns=s.sig_vars)
    floor_by_std = pd.DataFrame({k: floor.iloc[s.sets[k]].median() for k in core.STANDARD_TAGS if s.sets[k].size})
    cons = s.standard_consistency()
    cons_flat = cons.copy(); cons_flat.columns = [f'{a}_{b}' for a, b in cons.columns]
    io_xlsx.write_row_table(os.path.join(data_dir, 'StandardConsistency.xlsx'), cons_flat, 'Sheet1')
    dh = s.dh_slopes if s.dh_slopes is not None else s.downhole_slopes()
    io_xlsx.write_row_table(os.path.join(data_dir, 'DownholeSlopes.xlsx'), dh, 'Sheet1')
    log('StandardConsistency.xlsx and DownholeSlopes.xlsx written to data/')

    # ---- report
    report = os.path.join(out_dir, 'REPORT.html')
    write_report(report, s, kind, src, winfo, drift, excl, pd.DataFrame(cal_rows), pd.DataFrame(rec_rows),
                 floor_by_std, cons, dh)
    log(f'report: {report}')
    return report


def bundled_standards_path() -> str:
    base = getattr(sys, '_MEIPASS', None) or os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    for cand in (os.path.join(base, 'LAICPMS_Standard_Values.xlsx'),
                 os.path.join(os.path.dirname(base), 'LAICPMS_Standard_Values.xlsx'),
                 os.path.join(os.path.dirname(base), 'OG CODE', 'LAICPMS_Standard_Values.xlsx')):
        if os.path.isfile(cand):
            return cand
    raise FileNotFoundError('LAICPMS_Standard_Values.xlsx not found next to the program')


# ---------------------------------------------------------------- report
def write_report(path, s: ReductionSession, kind, src, winfo, drift, excl, cal_df, rec_df,
                 floor_by_std=None, cons=None, dh=None):
    e = html.escape
    zero_rows = [n for n, row in zip(s.display_names, s.final_result) if np.all(row == 0)]
    piv_b = rec_df.pivot(index='Analyte', columns='Standard', values='Bias%')
    piv_r = rec_df.pivot(index='Analyte', columns='Standard', values='RSD%_between')
    stds = [c for c in ('BHV', 'BCR', 'BIR', 'GSD', 'GSE', 'STH') if c in piv_b.columns]
    cal_ix = cal_df.set_index('Analyte')

    def fmt(v, nd=1):
        return '' if v is None or (isinstance(v, float) and not np.isfinite(v)) else f'{v:.{nd}f}'

    def flag(bias, rsd):
        cls = ''
        if (np.isfinite(bias) and abs(bias) > 10) or (np.isfinite(rsd) and rsd > 10):
            cls = ' class="bad"'
        elif (np.isfinite(bias) and abs(bias) > 5) or (np.isfinite(rsd) and rsd > 5):
            cls = ' class="warn"'
        return cls

    rows_html = []
    for a in cal_df.Analyte:
        r2 = float(cal_ix.loc[a, 'R2'])
        r2cls = ' class="bad"' if r2 < 0.99 else ''
        cells = [f'<td>{e(a)}</td><td{r2cls}>{r2:.4f}</td>']
        for st in stds:
            b = piv_b.loc[a, st] if a in piv_b.index else np.nan
            r = piv_r.loc[a, st] if a in piv_r.index else np.nan
            fl = float(floor_by_std.loc[a, st]) if (floor_by_std is not None and st in floor_by_std.columns and a in floor_by_std.index) else np.nan
            hide_bias = np.isfinite(fl) and fl > BIAS_FLOOR_PCT
            btxt = '<span class="muted">n/a</span>' if hide_bias else fmt(b)
            cells.append(f'<td{flag(np.nan if hide_bias else b, r)}>{btxt} / {fmt(r)} <span class="muted">({fmt(fl)})</span></td>')
        ex = cal_ix.loc[a, 'Excluded']
        cells.append(f'<td>{e(str(ex)) if ex != "None" else ""}</td>')
        rows_html.append('<tr>' + ''.join(cells) + '</tr>')

    cons_html = ''
    prov = [str(n) for n in s.std_vals[core.std_name_column(s.std_vals)].astype(str) if 'provisional' in str(n).lower()] if s.std_vals is not None else []
    prov_note = (f' <b>Provisional reference rows ({e(", ".join(prov))}) were derived from a previous run with this tool; '
                 f'their agreement is self-referential.</b>') if prov else ''
    if cons is not None and len(cons):
        keys = sorted({k for k, _ in cons.columns})
        keys = [k for k in ('BHV', 'BCR', 'BIR', 'GSD', 'GSE', 'STH', 'VE32', 'GOR_128') if k in keys]
        crow = []
        for a in cons.index:
            cells = [f'<td>{e(a)}</td>']
            for k in keys:
                v = cons.loc[a, (k, 'ratio')] if (k, 'ratio') in cons.columns else np.nan
                er = cons.loc[a, (k, 'err%')] if (k, 'err%') in cons.columns else np.nan
                cls = ''
                if np.isfinite(v) and np.isfinite(er) and er < BIAS_FLOOR_PCT:
                    dev = abs(100 * (v - 1))
                    cls = ' class="bad"' if dev > 2 * CONSISTENCY_FLAG_PCT else (' class="warn"' if dev > CONSISTENCY_FLAG_PCT else '')
                cells.append(f'<td{cls}>{fmt(v, 3)} <span class="muted">({fmt(er)})</span></td>' if np.isfinite(v) else '<td></td>')
            crow.append('<tr>' + ''.join(cells) + '</tr>')
        cons_html = (f'<h2>Consistency of the standards</h2><p>Apparent sensitivity of each standard relative to the '
                     f'calibration (1.000 = agrees with the calibrants), with its counting error % in brackets. '
                     f'A standard that is off for one element while the others agree points at that reference value '
                     f'or at a matrix effect; a standard that is off for every element points at its internal-standard '
                     f'oxide or its normaliser signal. Yellow &gt; {CONSISTENCY_FLAG_PCT:g} %, red &gt; {2 * CONSISTENCY_FLAG_PCT:g} % '
                     f'(only where the counting error is &lt; {BIAS_FLOOR_PCT:g} %). Table: <code>data/StandardConsistency.xlsx</code>.'
                     f'{prov_note}</p>'
                     f'<table><tr><th>Analyte</th>{"".join(f"<th>{k}</th>" for k in keys)}</tr>{"".join(crow)}</table>')

    dh_html = ''
    if dh is not None and len(dh):
        drow = ''.join(f'<tr><td>{e(a)}</td><td>{fmt(dh.loc[a, "slope_pct_per_10s"], 2)}</td>'
                       f'<td>{fmt(dh.loc[a, "spread_pct_per_10s"], 2)}</td></tr>' for a in dh.index)
        dh_html = (f'<h2>Down-hole fractionation (standards)</h2><p>Slope of analyte/normaliser inside the signal '
                   f'window, % per 10 s (median over the standard analyses; spread = robust sigma between analyses). '
                   f'An analysis integrated over a window shorter than the standards\' by &Delta;L has its ratios multiplied by '
                   f'exp(slope &times; &Delta;L / 2), which refers them to the standards\' window midpoint (on this run\'s '
                   f'standards that halves the 10-s-window bias to ~0.6 %). Table: <code>data/DownholeSlopes.xlsx</code>.</p>'
                   f'<table><tr><th>Analyte</th><th>slope %/10 s</th><th>spread</th></tr>{drow}</table>')

    short_html = ''
    if winfo.get('short_signal'):
        items = ', '.join(f'{e(s.display_names[i])} ({l:.0f} s)' for i, l in zip(winfo['short_signal'], winfo['short_signal_len']))
        dh = s.dh_factor
        corr_txt = ''
        if dh is not None and not np.allclose(dh, 1.0):
            mx = max(abs(100 * (dh[i] - 1)).max() for i in winfo['short_signal'])
            corr_txt = (f' Their ratios were referred to the standards\' window with the standards\' down-hole slopes '
                        f'(largest correction {mx:.1f} %; factors in <code>ReducedDataExport.xlsx / Downhole_Factor</code>).')
        short_html = (f'<li><b>{len(winfo["short_signal"])} analyses ablated for less than the standards</b> and were integrated '
                      f'over their own shorter window: {items}.{corr_txt}</li>')

    sets_txt = ', '.join(f'{k} ×{s.sets[k].size}' for k in core.STANDARD_TAGS if s.sets[k].size)
    n_auto = sum(len(el.auto_excluded) for el in s.cal_elements)
    by_norm = {}
    for el in s.cal_elements:
        by_norm.setdefault(s.norm_map[el.norm_key]['iso'], []).append(el.analyte)
    norm_html = '; '.join(f'<b>{e(k)}</b>: {e(" ".join(v))}' for k, v in by_norm.items())
    special = [(el.analyte, ' '.join(el.calibrants)) for el in s.cal_elements if list(el.calibrants) != list(core.PRIMARY_STDS)]
    special_html = (' Calibrants overridden for ' + ', '.join(f'{e(a)} ({e(c)})' for a, c in special) + '.') if special else ''
    doc = f"""<!DOCTYPE html><html><head><meta charset="utf-8"><title>Data check report</title>
<style>body{{font-family:Segoe UI,Arial,sans-serif;margin:24px;color:#1f2328}}table{{border-collapse:collapse;font-size:13px}}
td,th{{border:1px solid #d1d5db;padding:3px 8px;text-align:right}}th{{background:#f3f4f6}}td:first-child{{text-align:left}}
.warn{{background:#fef3c7}}.bad{{background:#fecaca}}.muted{{color:#6b7280;font-size:11px}}img{{max-width:100%}}code{{background:#f3f4f6;padding:1px 4px}}</style></head><body>
<h1>LA-ICP-MS data check</h1>
<p>{e(datetime.now().strftime('%Y-%m-%d %H:%M'))} &nbsp; input: <code>{e(src)}</code> ({kind})</p>
<h2>Run</h2><ul>
<li>{s.n_samples} analyses, {len(s.sig_vars)} analytes: {e(' '.join(s.sig_vars))}</li>
<li>standards found: {e(sets_txt)}</li>
<li>dwell per isotope: {e(', '.join(f'{v} {d * 1000:g} ms' for v, d in zip(s.sig_vars, s.dwell_s)) if s.dwell_s is not None else 'unknown')}</li></ul>
<h2>Choices made</h2><ul>
<li>Background: 1 s to laser-on − 2 s. Signal: laser-on + {SIGNAL_DELAY_S:g} s for {winfo['signal_length_s']:.1f} s (laser-on median {winfo['laser_on_median_s']:.1f} s, standards' ablation median {winfo['ablation_median_s']:.1f} s). Second background: {winfo['second_background']}.</li>
{short_html}
<li>Drift: {e(drift[2])}.</li>
<li>Calibration: BHV, BCR, BIR vs GeoRem values. Normaliser chosen per analyte from the standards (lowest replicate RSD and cross-standard spread in quadrature, <code>data/NormaliserChoice.xlsx</code>): {norm_html}; slope = inverse-variance weighted mean of log(target/measured) over the calibrant replicates; {n_auto} replicate values dropped for counting error &gt; {core.MAX_CAL_ERR_PCT:g} %.{special_html}</li>
<li>Auto-excluded calibrant replicates (&gt; 4 robust σ and &gt; 10 % from their own standard's other replicates): {e('; '.join(f'{k}: {" ".join(v)}' for k, v in excl.items()) if excl else 'none')}.</li>
<li>Internal-standard oxide wt% were autofilled for BHV, BCR, BIR, GSD, GSE, StHS, VE32 and GOR-128 only. <b>{len(zero_rows)} samples have no oxide value and therefore read 0</b>: fill <code>data/Intervals.xlsx</code> and re-run the reduction.</li>
</ul>
<h2>Calibration quality</h2>
<p>Per analyte: R², then for each standard <i>bias % / between-replicate RSD % (counting-statistics floor %)</i>. The floor is the median counting error of that standard's replicates for the analyte/normaliser ratio; an RSD near the floor cannot be improved by the reduction. Bias is not shown where the floor exceeds {BIAS_FLOOR_PCT:g} %. Yellow &gt; 5 %, red &gt; 10 % (or R² &lt; 0.99).</p>
<table><tr><th>Analyte</th><th>R²</th>{''.join(f'<th>{st}</th>' for st in stds)}<th>excluded</th></tr>
{''.join(rows_html)}</table>
{cons_html}
{dh_html}
<h2>Plots</h2>
<p><a href="plots/calibration_curves/index.png">Calibration curves</a> (one PNG per analyte in <code>plots/calibration_curves/</code>)</p>
<img src="plots/calibration_curves/index.png">
<p><a href="plots/standard_recovery/index.png">Standard recovery</a> (one PNG per analyte in <code>plots/standard_recovery/</code>; bars = ±2 SEM per analysis)</p>
<img src="plots/standard_recovery/index.png">
<h2>Files</h2><ul>
<li><code>data/ReducedDataExport.xlsx</code> – FinalResult, Signal_SEM_Norm, Counting_Err_Pct, snapshots of all settings</li>
<li><code>data/DriftCorrected_All.xlsx</code> – drift-corrected background-subtracted signals</li>
<li><code>data/StandardConsistency.xlsx</code>, <code>data/DownholeSlopes.xlsx</code> – the two tables above</li>
<li><code>data/Intervals.xlsx</code>, <code>DriftSelections.xlsx</code>, <code>CalibrationSelections.xlsx</code> – the choices above, editable; re-run with "Re-run a *_checked folder" after editing</li>
</ul></body></html>"""
    with open(path, 'w', encoding='utf-8') as f:
        f.write(doc)
