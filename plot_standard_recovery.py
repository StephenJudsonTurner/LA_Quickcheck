"""Standard recovery plots: reference value vs. reduced value, linear axes, for
BHV / BCR / BIR / GSD on every analyte, with error bars from the internal
uncertainty of each individual analysis.

    python plot_standard_recovery.py <data file.xlsx> [--std FILE] [--out DIR] [--tag NAME]
                                     [--errors sem|rsd] [--k 2]

x axis : reduced value from FinalResult (the same numbers ReducedDataExport.xlsx holds)
y axis : reference value from LAICPMS_Standard_Values.xlsx (BHV/BCR/BIR from the
         GeoRem or Harvard row the analyte's calibration used; GSD from the GSD1G row)
x error: +/- k * (internal uncertainty %) * reduced value, where the internal
         uncertainty is the tool's per-analysis SEM% of the analyte/normaliser ratio
         over the central 50% of the signal window (--errors sem, default, k=2 ->
         2 SEM) or the RSD% of the same frames (--errors rsd).
Also drawn: the 1:1 line; points excluded from the calibration fit are marked
with a red X.  A stats box lists, per standard, n, mean, between-replicate
%RSD, the reference value and %bias, i.e. the tool's secondary-standards table.
`recovery_stats.xlsx` tabulates the same numbers for all analytes.

The source data folder is not modified (inputs are copied to a temp folder).
"""
from __future__ import annotations

import argparse
import os
import shutil
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lasercal.session import ReductionSession, UIState, IntervalsCreated  # noqa: E402
from lasercal import core, io_xlsx  # noqa: E402
from plot_calibrations import stage_copy, contact_sheet, INK, MUTED, EXCL, SURFACE  # noqa: E402

# fixed order, one hue + one marker per standard (identity is carried by marker + label too)
STANDARDS = [  # key in session.sets, label, standards-table match, marker, colour
    ('BHV', 'BHV', 'BHV', 'o', '#2563eb'),
    ('BCR', 'BCR', 'BCR', '^', '#d97706'),
    ('BIR', 'BIR', 'BIR', 's', '#059669'),
    ('GSD', 'GSD', 'GSD1G', 'D', '#7c3aed'),
]


def reference_value(std_vals: pd.DataFrame, analyte: str, match: str, set_name: str) -> float:
    """Row lookup as the MATLAB summary table does it: primaries need the set tag."""
    if analyte not in std_vals.columns:
        return np.nan
    names = [core._norm_str(s) for s in std_vals[core.std_name_column(std_vals)].astype(str)]
    key = core._norm_str(match)
    if match in ('BHV', 'BCR', 'BIR'):
        hits = [k for k, n in enumerate(names) if key in n and core._norm_str(set_name) in n]
    else:
        hits = [k for k, n in enumerate(names) if n == key or key in n]
    if not hits:
        return np.nan
    v = pd.to_numeric(pd.Series([std_vals[analyte].iloc[hits[0]]]), errors='coerce').iloc[0]
    return float(v) if pd.notna(v) else np.nan


def plot_one(s, el, j, err_pct, k, err_name, out_png):
    final = s.final_result[:, j]
    fig, ax = plt.subplots(figsize=(8, 8), dpi=130)
    fig.patch.set_facecolor(SURFACE); ax.set_facecolor(SURFACE)
    ax.grid(True, color='#e5e7eb', lw=0.8)
    for sp in ax.spines.values():
        sp.set_color('#d1d5db')

    xs, ys, rows = [], [], []
    for key, label, match, marker, colour in STANDARDS:
        idx = s.sets.get(key, np.array([], int))
        if idx.size == 0:
            continue
        ref = reference_value(s.std_vals, el.analyte, match, el.set_name)
        vals = final[idx]
        errs = k * np.abs(err_pct[idx, j]) / 100.0 * np.abs(vals)
        names = [s.display_names[i] for i in idx]
        exc = np.array([n in el.excluded for n in names], bool)
        fin = np.isfinite(vals)
        if np.isfinite(ref) and fin.any():
            ax.errorbar(vals[fin], np.full(fin.sum(), ref), xerr=np.where(np.isfinite(errs[fin]), errs[fin], 0),
                        fmt='none', ecolor=colour, elinewidth=1, capsize=3, alpha=0.8, zorder=2)
            ax.scatter(vals[fin], np.full(fin.sum(), ref), marker=marker, s=55, facecolor='none',
                       edgecolor=colour, linewidth=1.6, zorder=3, label=f'{label} (n={int(fin.sum())})')
            if (exc & fin).any():
                ax.scatter(vals[exc & fin], np.full((exc & fin).sum(), ref), marker='x', s=90, color=EXCL,
                           linewidth=1.8, zorder=4)
            ax.annotate(f'{label}\nref = {ref:.4g}', xy=(np.nanmax(vals[fin]), ref), xytext=(10, 0),
                        textcoords='offset points', fontsize=8, va='center', ha='left', color=INK, zorder=5)
            xs.extend(vals[fin].tolist()); ys.append(ref)
        mu, sd = core.nanmean(vals), core.nanstd(vals)
        n = int(fin.sum())
        rsd = 100 * sd / mu if (np.isfinite(mu) and mu != 0 and n > 1) else np.nan
        bias = 100 * (mu - ref) / ref if (np.isfinite(mu) and np.isfinite(ref) and ref != 0) else np.nan
        rows.append({'Analyte': el.analyte, 'Standard': label, 'N': n, 'Mean': mu, 'RSD%_between': rsd,
                     'Reference': ref, 'Mean/Ref': mu / ref if np.isfinite(bias) else np.nan, 'Bias%': bias,
                     'MeanInternal%': core.nanmean(err_pct[idx, j]), 'nExcludedFromFit': int(exc.sum())})

    if xs and ys:
        lo = min(min(xs), min(ys)); hi = max(max(xs), max(ys))
        pad = 0.08 * (hi - lo if hi > lo else abs(hi) or 1.0)
        lo, hi = lo - pad, hi + pad
        ax.plot([lo, hi], [lo, hi], ls='--', color='#9ca3af', lw=1.2, zorder=1, label='1:1')
        ax.set_xlim(lo, hi + 0.25 * (hi - lo)); ax.set_ylim(lo, hi)

    ax.set_xlabel(f'Reduced value  {el.analyte}   (error bars: ±{k:g} {err_name} of each analysis)', color=INK)
    ax.set_ylabel(f'Reference value  {el.analyte}   (standards table; primaries: {el.set_name})', color=INK)
    ax.set_title(f'Standard recovery: {el.analyte}   [normalised to {el.norm_key}]', color=INK, fontsize=12)

    lines = [f'{"Std":<4} {"n":>2} {"mean":>9} {"RSD%":>6} {"ref":>9} {"bias%":>6}']
    for r in rows:
        f = lambda v, w: f'{v:>{w}.4g}' if np.isfinite(v) else f'{"n/a":>{w}}'
        lines.append(f'{r["Standard"]:<4} {r["N"]:>2} {f(r["Mean"], 9)} {f(r["RSD%_between"], 6)} '
                     f'{f(r["Reference"], 9)} {f(r["Bias%"], 6)}')
    if el.excluded:
        lines.append('excluded from fit: ' + ' '.join(el.excluded))
    ax.text(0.02, 0.98, '\n'.join(lines), transform=ax.transAxes, va='top', ha='left', fontsize=8,
            family='monospace', color=INK,
            bbox=dict(boxstyle='round,pad=0.4', facecolor='white', edgecolor='#d1d5db', alpha=0.92), zorder=6)
    ax.legend(loc='lower right', fontsize=8, frameon=True, edgecolor='#d1d5db')
    fig.tight_layout()
    fig.savefig(out_png, facecolor=SURFACE)
    plt.close(fig)
    return rows


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('data_file')
    ap.add_argument('--std', help='standards workbook (default: LAICPMS_Standard_Values.xlsx next to this script)')
    ap.add_argument('--out'); ap.add_argument('--tag')
    ap.add_argument('--errors', choices=['sem', 'rsd'], default='sem',
                    help='per-analysis uncertainty: SEM%% of the plateau ratio (default) or RSD%%')
    ap.add_argument('--k', type=float, default=2.0, help='multiplier for the error bars (default 2 -> 2 SEM)')
    a = ap.parse_args(argv)

    here = os.path.dirname(os.path.abspath(__file__))
    tag = a.tag or os.path.basename(os.path.dirname(os.path.abspath(a.data_file))) + '_recovery'
    out = a.out or os.path.join(here, 'plots', tag)
    os.makedirs(out, exist_ok=True)

    staged = stage_copy(a.data_file)
    s = ReductionSession(staged, std_file=a.std, ui=UIState())
    try:
        s.load()
    except IntervalsCreated as e:
        print(e); return 2
    s.run_drift(); s.calibrate()
    err_pct = core.compute_sem_backcorr_normalized(s.windows, s.intervals, s.sig_vars, s.norm_map,
                                                   s.load_or_init_cal_prefs(), stat=a.errors)
    err_name = 'SEM' if a.errors == 'sem' else 'RSD'

    pngs, rows = [], []
    for j, el in enumerate(s.cal_elements):
        png = os.path.join(out, f'{el.analyte}.png')
        rows += plot_one(s, el, j, err_pct, a.k, err_name, png)
        pngs.append(png)
    io_xlsx.write_table(os.path.join(out, 'recovery_stats.xlsx'), pd.DataFrame(rows), 'Sheet1')
    contact_sheet(pngs, os.path.join(out, 'index.png'))
    with open(os.path.join(out, 'SOURCE.txt'), 'w', encoding='utf-8') as f:
        f.write(f'data file: {os.path.abspath(a.data_file)}\nstandards: {s.std_file}\n'
                f'error bars: {a.k:g} x {err_name}% of each analysis (plateau ratio statistic)\n'
                f'samples: {s.n_samples}, analytes: {len(s.sig_vars)}\n')
    shutil.rmtree(os.path.dirname(staged), ignore_errors=True)
    print(f'wrote {len(pngs)} PNGs + index.png + recovery_stats.xlsx to {out}')
    return 0


if __name__ == '__main__':
    sys.exit(main())
