"""Automatic first-pass reduction ("data checker").

Given either a folder of Thermo Qtegra per-analysis CSVs or a MATLAB-style
data folder (data workbook + RunOrder.xlsx), choose reduction parameters
from the data itself, run the reduction, and write the standard plots and an
HTML report.  Every choice is written to the output folder in the same XLSX
files the tool normally reads, so a second pass can edit them by hand.

Rules used for the automatic choices (see choose_* functions):

* windows       laser-on / laser-off from the 25Mg trace of each analysis;
                background = [1 s, laser-on - 2 s]; signal starts 1 s after
                laser-on and lasts min(20 s, ablation length - 5 s), the same
                length for every sample; a second background only if >= 8 s
                of clean gas blank follow the ablation.
* drift         standard = BHV if present, else the most frequent of BCR/GSD/BIR.
                Model from how that standard was run: interpolation (pchip)
                when it was run in consecutive pairs or more, quadratic
                (poly2) when run singly with >= 4 brackets, linear with 2-3
                brackets, none with 1.
* calibration   normalise to Al (Ca for 27Al) when Al was measured, else Ca;
                intercept forced through zero; GeoRem values.  A calibrant
                replicate that departs from its own standard's other
                replicates by > 4 robust sigma and > 10 % is excluded
                (never more than a third of a standard's points) and listed.
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

SIGNAL_DELAY_S = 1.0   # signal starts this long after laser-on
SIGNAL_CAP_S = 20.0    # and lasts min(this, ablation length - 5 s)


def _log_default(msg: str):
    print(msg, flush=True)


# ---------------------------------------------------------------- input detection
def detect_input(path: str) -> Tuple[str, str]:
    """-> ('qtegra', folder) | ('xlsx', data workbook path)."""
    if os.path.isdir(path):
        csvs = qtegra.list_files(path)
        if csvs:
            return 'qtegra', path
        xl = [f for f in os.listdir(path) if f.lower().endswith('.xlsx') and f not in
              ('RunOrder.xlsx', 'Intervals.xlsx', 'DriftSelections.xlsx', 'CalibrationSelections.xlsx',
               'SampleWindows.xlsx', 'DriftCorrected_All.xlsx', 'ReducedDataExport.xlsx', 'SourceFiles.xlsx',
               'LAICPMS_Standard_Values.xlsx')]
        if len(xl) == 1 and os.path.isfile(os.path.join(path, 'RunOrder.xlsx')):
            return 'xlsx', os.path.join(path, xl[0])
        raise ValueError('Folder must contain Qtegra CSV files, or exactly one data workbook plus RunOrder.xlsx.')
    if path.lower().endswith('.xlsx'):
        return 'xlsx', path
    raise ValueError(f'Unrecognised input: {path}')


# ---------------------------------------------------------------- choices
def choose_windows(windows: List[pd.DataFrame], thresh: float = core.START_THRESH,
                   iso: str = 'x25Mg') -> Tuple[pd.DataFrame, Dict]:
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
    # Window sweep on the 09_21_26_50um run (PythonTool/sweep_signal_window.py): precision
    # flattens beyond ~20 s and the delay after laser-on makes no difference, so:
    sig_len = float(np.nanmin(np.minimum(SIGNAL_CAP_S, ablation - 5.0)))
    sig_len = max(sig_len, 4.0)
    tail = end_t - off_t
    use_back2 = bool(np.nanmin(tail) >= 12.0)
    iv = pd.DataFrame(np.nan, index=range(len(windows)), columns=core.INTERVAL_COLUMNS)
    iv['back_start1'] = 1.0
    iv['back_stop1'] = on_t - 2.0
    iv['signal_start'] = on_t + SIGNAL_DELAY_S
    iv['signal_stop'] = on_t + SIGNAL_DELAY_S + sig_len
    if use_back2:
        iv['back_start2'] = off_t + 4.0
        iv['back_stop2'] = end_t - 0.5
    info = {'laser_on_median_s': float(np.nanmedian(on_t)), 'ablation_median_s': float(np.nanmedian(ablation)),
            'signal_length_s': sig_len, 'second_background': use_back2,
            'n_no_laser_on': int(np.isnan(on_t).sum())}
    return iv, info


def choose_drift(sets: Dict[str, np.ndarray]) -> Tuple[str, str, str]:
    """-> (standard, method, reason)."""
    order = ['BHV', 'BCR', 'GSD', 'BIR']
    present = [k for k in order if sets.get(k, np.array([])).size > 0]
    if not present:
        return 'BHV', 'none', 'no drift standard found; no drift correction'
    std = 'BHV' if 'BHV' in present else max(present, key=lambda k: sets[k].size)
    groups = core.group_consecutive(sets[std])
    n_groups = len(groups)
    paired = sum(1 for g in groups if g.size >= 2) >= max(2, n_groups // 2)
    if n_groups < 2:
        return std, 'none', f'{std} run only once; no drift correction possible'
    if paired:
        return std, 'pchip', f'{std} run in consecutive pairs ({n_groups} brackets): interpolate through bracket means'
    if n_groups >= 4:
        return std, 'poly2', f'{std} run singly in {n_groups} brackets: quadratic fit (an interpolant would pin every {std} to its mean)'
    return std, 'linear', f'{std} run singly in {n_groups} brackets: linear'


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
        rows.append({'NormElement': key, 'ForceInterceptZero': True, 'StandardSet': 'GeoRem', 'CalExclusions': 'None'})
    return pd.DataFrame(rows, index=sig_vars, dtype=object)


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
        ', '.join(f'{k} x{sets[k].size}' for k in ('BHV', 'BCR', 'BIR', 'GSD', 'STH', 'VE32', 'GOR_128') if sets[k].size))

    iv, winfo = choose_windows(s.windows)
    iv.index = s.display_names
    io_xlsx.write_intervals(s.intervals_path, iv)
    log(f'windows: laser-on at {winfo["laser_on_median_s"]:.1f} s, ablation {winfo["ablation_median_s"]:.1f} s; '
        f'signal = laser-on + {SIGNAL_DELAY_S:g} s for {winfo["signal_length_s"]:.1f} s; second background: {winfo["second_background"]}')

    dstd, dmeth, dreason = choose_drift(sets)
    io_xlsx.write_drift_selections(s.drift_prefs_path, pd.DataFrame(
        {'Standard': dstd, 'Interpolation': dmeth, 'Exclusions': ''}, index=s.sig_vars, dtype=object))
    log(f'drift: {dreason}')
    cal = choose_calibration(s.sig_vars, s.norm_map)
    io_xlsx.write_calibration_selections(s.cal_prefs_path, cal)
    log(f'calibration: normalised to {cal.NormElement.mode()[0]} (Ca for the Al isotope), intercept 0, GeoRem')

    # ---- reduce, find outliers, reduce again
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
    skip = {'RunOrder.xlsx', 'Intervals.xlsx', 'DriftSelections.xlsx', 'CalibrationSelections.xlsx', 'SampleWindows.xlsx',
            'DriftCorrected_All.xlsx', 'ReducedDataExport.xlsx', 'SourceFiles.xlsx', 'LAICPMS_Standard_Values.xlsx'}
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
    winfo = {'laser_on_median_s': float(np.nan), 'ablation_median_s': float(np.nan),
             'signal_length_s': float(np.nanmedian(iv['signal_stop'] - iv['signal_start'])),
             'second_background': bool(iv['back_start2'].notna().any())}
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

    # ---- report
    report = os.path.join(out_dir, 'REPORT.html')
    write_report(report, s, kind, src, winfo, drift, excl, pd.DataFrame(cal_rows), pd.DataFrame(rec_rows))
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
def write_report(path, s: ReductionSession, kind, src, winfo, drift, excl, cal_df, rec_df):
    e = html.escape
    zero_rows = [n for n, row in zip(s.display_names, s.final_result) if np.all(row == 0)]
    piv_b = rec_df.pivot(index='Analyte', columns='Standard', values='Bias%')
    piv_r = rec_df.pivot(index='Analyte', columns='Standard', values='RSD%_between')
    stds = [c for c in ('BHV', 'BCR', 'BIR', 'GSD') if c in piv_b.columns]

    def fmt(v):
        return '' if v is None or (isinstance(v, float) and not np.isfinite(v)) else f'{v:.1f}'

    def flag(bias, rsd):
        cls = ''
        if np.isfinite(bias) and abs(bias) > 10 or (np.isfinite(rsd) and rsd > 10):
            cls = ' class="bad"'
        elif np.isfinite(bias) and abs(bias) > 5 or (np.isfinite(rsd) and rsd > 5):
            cls = ' class="warn"'
        return cls

    rows_html = []
    for a in cal_df.Analyte:
        r2 = float(cal_df.set_index('Analyte').loc[a, 'R2'])
        r2cls = ' class="bad"' if r2 < 0.99 else ''
        cells = [f'<td>{e(a)}</td><td{r2cls}>{r2:.4f}</td>']
        for st in stds:
            b = piv_b.loc[a, st] if a in piv_b.index else np.nan
            r = piv_r.loc[a, st] if a in piv_r.index else np.nan
            cells.append(f'<td{flag(b, r)}>{fmt(b)} / {fmt(r)}</td>')
        ex = cal_df.set_index('Analyte').loc[a, 'Excluded']
        cells.append(f'<td>{e(str(ex)) if ex != "None" else ""}</td>')
        rows_html.append('<tr>' + ''.join(cells) + '</tr>')

    sets_txt = ', '.join(f'{k} ×{s.sets[k].size}' for k in ('BHV', 'BCR', 'BIR', 'GSD', 'STH', 'VE32', 'GOR_128') if s.sets[k].size)
    doc = f"""<!DOCTYPE html><html><head><meta charset="utf-8"><title>Data check report</title>
<style>body{{font-family:Segoe UI,Arial,sans-serif;margin:24px;color:#1f2328}}table{{border-collapse:collapse;font-size:13px}}
td,th{{border:1px solid #d1d5db;padding:3px 8px;text-align:right}}th{{background:#f3f4f6}}td:first-child{{text-align:left}}
.warn{{background:#fef3c7}}.bad{{background:#fecaca}}img{{max-width:100%}}code{{background:#f3f4f6;padding:1px 4px}}</style></head><body>
<h1>LA-ICP-MS data check</h1>
<p>{e(datetime.now().strftime('%Y-%m-%d %H:%M'))} &nbsp; input: <code>{e(src)}</code> ({kind})</p>
<h2>Run</h2><ul>
<li>{s.n_samples} analyses, {len(s.sig_vars)} analytes: {e(' '.join(s.sig_vars))}</li>
<li>standards found: {e(sets_txt)}</li></ul>
<h2>Choices made</h2><ul>
<li>Background: 1 s to laser-on − 2 s. Signal: laser-on + 1 s for {winfo['signal_length_s']:.1f} s (laser-on median {winfo['laser_on_median_s']:.1f} s, ablation median {winfo['ablation_median_s']:.1f} s). Second background: {winfo['second_background']}.</li>
<li>Drift: {e(drift[2])}.</li>
<li>Calibration: BHV, BCR, BIR vs GeoRem values, normalised to {e(s.norm_map['Al']['iso'] if s.norm_map['Al']['iso'] in s.sig_vars else s.norm_map['Ca']['iso'])} (27Al to Ca), intercept forced to zero.</li>
<li>Auto-excluded calibrant replicates (&gt; 4 robust σ and &gt; 10 % from their own standard's other replicates): {e('; '.join(f'{k}: {" ".join(v)}' for k, v in excl.items()) if excl else 'none')}.</li>
<li>Internal-standard oxide wt% were autofilled for BHV, BCR, BIR, GSD, StHS, VE32 and GOR-128 only. <b>{len(zero_rows)} samples have no oxide value and therefore read 0</b>: fill <code>data/Intervals.xlsx</code> and re-run the reduction.</li>
</ul>
<h2>Calibration quality</h2>
<p>Per analyte: R², then for each standard <i>bias % / between-replicate RSD %</i>. Yellow &gt; 5 %, red &gt; 10 % (or R² &lt; 0.99).</p>
<table><tr><th>Analyte</th><th>R²</th>{''.join(f'<th>{st}</th>' for st in stds)}<th>excluded</th></tr>
{''.join(rows_html)}</table>
<h2>Plots</h2>
<p><a href="plots/calibration_curves/index.png">Calibration curves</a> (one PNG per analyte in <code>plots/calibration_curves/</code>)</p>
<img src="plots/calibration_curves/index.png">
<p><a href="plots/standard_recovery/index.png">Standard recovery</a> (one PNG per analyte in <code>plots/standard_recovery/</code>; bars = ±2 SEM per analysis)</p>
<img src="plots/standard_recovery/index.png">
<h2>Files</h2><ul>
<li><code>data/ReducedDataExport.xlsx</code> – FinalResult, Signal_SEM_Norm, snapshots of all settings</li>
<li><code>data/DriftCorrected_All.xlsx</code> – drift-corrected background-subtracted signals</li>
<li><code>data/Intervals.xlsx</code>, <code>DriftSelections.xlsx</code>, <code>CalibrationSelections.xlsx</code> – the choices above, editable; re-run with <code>run_reduction</code> after editing</li>
</ul></body></html>"""
    with open(path, 'w', encoding='utf-8') as f:
        f.write(doc)
