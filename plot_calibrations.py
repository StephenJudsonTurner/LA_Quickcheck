"""Primary-standard calibration curves (BHV / BCR / BIR) for every analyte, as PNG.

    python plot_calibrations.py <data file.xlsx> [--std FILE] [--out DIR] [--tag NAME]

One PNG per analyte reproducing the Calibration tab of the MATLAB tool:
log-log scatter of measured (drift-corrected, normalised) vs. standard target
(normalised the same way), every point labelled with its sample name, excluded
points drawn as red X, the fitted line from the included points, and a stats
box (slope, intercept, R^2, n, normaliser, standard set, exclusions).
A `calibration_stats.xlsx` table with the same numbers is written next to
the PNGs, plus an `index.png` contact sheet of all curves.

The source data folder is never modified: its input files are copied to a
temp folder and the pipeline runs there (the MATLAB tool re-saves
Intervals.xlsx on launch; this avoids doing that to your originals).
"""
from __future__ import annotations

import argparse
import os
import shutil
import sys
import tempfile

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lasercal.session import ReductionSession, UIState, IntervalsCreated  # noqa: E402
from lasercal import io_xlsx  # noqa: E402

INPUT_FILES = ('RunOrder.xlsx', 'Intervals.xlsx', 'DriftSelections.xlsx', 'CalibrationSelections.xlsx',
               'SampleWindows.xlsx')

INK = '#1f2328'          # text
MUTED = '#6b7280'        # secondary text / grid
POINT = '#2563eb'        # included primary-standard points (single series)
FIT = '#374151'          # fit line
EXCL = '#dc2626'         # excluded marker
SURFACE = '#ffffff'


def stage_copy(data_file: str) -> str:
    """Copy the data workbook and the input XLSX files to a temp folder."""
    src = os.path.dirname(os.path.abspath(data_file))
    work = tempfile.mkdtemp(prefix='lasercal_plots_')
    shutil.copy2(data_file, os.path.join(work, os.path.basename(data_file)))
    for f in INPUT_FILES:
        p = os.path.join(src, f)
        if os.path.isfile(p):
            shutil.copy2(p, os.path.join(work, f))
    return os.path.join(work, os.path.basename(data_file))


def std_group(name: str) -> str:
    low = name.lower()
    return 'BHV' if 'bhv' in low else 'BCR' if 'bcr' in low else 'BIR' if 'bir' in low else '?'


def plot_one(el, s, out_png: str):
    x, y = el.x_plot, el.y_plot
    names = el.names_plot
    excluded = np.array([n in el.excluded for n in names], bool)
    plottable = (x > 0) & (y > 0)
    n_hidden = int((~plottable).sum())

    fig, ax = plt.subplots(figsize=(8, 8), dpi=130)
    fig.patch.set_facecolor(SURFACE)
    ax.set_facecolor(SURFACE)
    ax.set_xscale('log'); ax.set_yscale('log')
    ax.grid(True, which='major', color='#e5e7eb', lw=0.8)
    ax.grid(True, which='minor', color='#f3f4f6', lw=0.5)
    for sp in ax.spines.values():
        sp.set_color('#d1d5db')

    inc = plottable & ~excluded
    exc = plottable & excluded
    ax.scatter(x[inc], y[inc], s=55, facecolor='none', edgecolor=POINT, linewidth=1.6, zorder=3,
               label=f'included (n={int(inc.sum())})')
    if exc.any():
        ax.scatter(x[exc], y[exc], s=55, facecolor='none', edgecolor=MUTED, linewidth=1.2, zorder=3)
        ax.scatter(x[exc], y[exc], s=90, marker='x', color=EXCL, linewidth=1.8, zorder=4,
                   label=f'excluded (n={int(exc.sum())})')

    # fit line over the plotted x-range (like updateCalibrationPlot)
    if np.isfinite(el.slope) and plottable.any():
        xx = np.geomspace(x[plottable].min(), x[plottable].max(), 100)
        yy = el.slope * xx + el.intercept
        ok = np.isfinite(yy) & (yy > 0)
        if ok.any():
            ax.plot(xx[ok], yy[ok], color=FIT, lw=2, zorder=2, label='fit')

    # axis limits with room on the right for the label columns
    if plottable.any():
        xlo, xhi = x[plottable].min(), x[plottable].max()
        ylo, yhi = y[plottable].min(), y[plottable].max()
        ax.set_xlim(xlo / 1.15, xhi * 1.6)
        ax.set_ylim(ylo / 1.35, yhi * 1.35)

    # one label per standard, beside its cluster (all points of a standard share the target y)
    for g in ('BHV', 'BCR', 'BIR'):
        m = np.flatnonzero(np.array([std_group(n) == g for n in names]) & plottable)
        if m.size == 0:
            continue
        n_exc = int(excluded[m].sum())
        txt = f'{g}  (n={m.size}' + (f', {n_exc} excluded' if n_exc else '') + f')\ntarget = {y[m][0]:.4g}'
        ax.annotate(txt, xy=(x[m].max(), y[m][0]), xytext=(10, 0), textcoords='offset points',
                    fontsize=8, va='center', ha='left', color=INK, zorder=5)

    norm_iso = s.norm_map[el.norm_key]['iso']
    ax.set_xlabel(f'Measured, drift-corrected  {el.analyte} / {norm_iso}', color=INK)
    ax.set_ylabel(f'Standard target  {el.analyte} / {norm_iso}  (set: {el.set_name})', color=INK)
    ax.set_title(f'Calibration: {el.analyte}   [{el.set_name}, normalised to {el.norm_key}]', color=INK, fontsize=12)

    lines = [f'Slope: {el.slope:.6g}',
             'Intercept: 0 (forced)' if el.force_zero else f'Intercept: {el.intercept:.6g}',
             f'R²: {el.r2:.5f}' if np.isfinite(el.r2) else 'R²: n/a',
             f'n included: {int((~excluded).sum())}   n excluded: {int(excluded.sum())}']
    if el.excluded:
        lines.append('Excluded: ' + ' '.join(el.excluded))
    if n_hidden:
        lines.append(f'{n_hidden} point(s) ≤ 0 not shown on log axes')
    ax.text(0.02, 0.98, '\n'.join(lines), transform=ax.transAxes, va='top', ha='left', fontsize=8.5,
            color=INK, bbox=dict(boxstyle='round,pad=0.4', facecolor='white', edgecolor='#d1d5db', alpha=0.9),
            zorder=6)
    ax.legend(loc='lower right', fontsize=8, frameon=True, edgecolor='#d1d5db')
    fig.tight_layout()
    fig.savefig(out_png, facecolor=SURFACE)
    plt.close(fig)


def contact_sheet(pngs, out_png):
    import matplotlib.image as mpimg
    n = len(pngs)
    cols = 6
    rows = int(np.ceil(n / cols))
    fig, axes = plt.subplots(rows, cols, figsize=(cols * 3.2, rows * 3.2), dpi=80)
    for ax, p in zip(axes.ravel(), pngs):
        ax.imshow(mpimg.imread(p)); ax.set_title(os.path.splitext(os.path.basename(p))[0], fontsize=8)
    for ax in axes.ravel():
        ax.axis('off')
    fig.tight_layout()
    fig.savefig(out_png)
    plt.close(fig)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('data_file')
    ap.add_argument('--std', help='standards workbook (default: LAICPMS_Standard_Values.xlsx next to this script)')
    ap.add_argument('--out', help='output folder (default: PythonTool/plots/<tag>)')
    ap.add_argument('--tag', help='sub-folder name under plots/ (default: data folder name)')
    a = ap.parse_args(argv)

    here = os.path.dirname(os.path.abspath(__file__))
    tag = a.tag or os.path.basename(os.path.dirname(os.path.abspath(a.data_file)))
    out = a.out or os.path.join(here, 'plots', tag)
    os.makedirs(out, exist_ok=True)

    staged = stage_copy(a.data_file)
    s = ReductionSession(staged, std_file=a.std, ui=UIState())
    try:
        s.load()
    except IntervalsCreated as e:
        print(e); return 2
    s.run_drift()
    s.calibrate()

    pngs, rows = [], []
    for el in s.cal_elements:
        png = os.path.join(out, f'{el.analyte}.png')
        plot_one(el, s, png)
        pngs.append(png)
        n_exc = len([n for n in el.names_plot if n in el.excluded])
        rows.append({'Analyte': el.analyte, 'Normaliser': s.norm_map[el.norm_key]['iso'], 'StandardSet': el.set_name,
                     'InterceptForcedZero': el.force_zero, 'Slope': el.slope, 'Intercept': el.intercept,
                     'R2': el.r2, 'nIncluded': len(el.names_plot) - n_exc, 'nExcluded': n_exc,
                     'Excluded': ' '.join(el.excluded) if el.excluded else 'None'})
    io_xlsx.write_table(os.path.join(out, 'calibration_stats.xlsx'), pd.DataFrame(rows), 'Sheet1')
    contact_sheet(pngs, os.path.join(out, 'index.png'))
    with open(os.path.join(out, 'SOURCE.txt'), 'w', encoding='utf-8') as f:
        f.write(f'data file: {os.path.abspath(a.data_file)}\nstandards: {s.std_file}\n'
                f'inputs: {", ".join(INPUT_FILES)} from the same folder, copied at run time\n'
                f'samples: {s.n_samples}, analytes: {len(s.sig_vars)}\n')
    shutil.rmtree(os.path.dirname(staged), ignore_errors=True)
    print(f'wrote {len(pngs)} PNGs + index.png + calibration_stats.xlsx to {out}')
    return 0


if __name__ == '__main__':
    sys.exit(main())
