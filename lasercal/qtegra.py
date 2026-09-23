"""Adapter for Thermo Qtegra per-analysis CSV exports (iCAP RQ).

Each analysis is one CSV: line 1 is ``<sample name>:<timestamp>;``, a dozen
instrument-metadata lines follow, then a ``Time,23Na,25Mg,...`` header, a
dwell-time line, and the time-resolved intensities.

``build_data_folder`` turns a folder of those files into the input set the
reduction tool expects:

* ``<name>_data.xlsx``  all analyses concatenated (Time = seconds since the
  first analysis started + time inside the file), MATLAB-style column names
* ``RunOrder.xlsx``      sample names in acquisition order
* ``SampleWindows.xlsx`` explicit window per sample (first/last row, laser-on
  row) so the tool uses each file as its own window instead of the
  threshold search, which cannot separate files that end during washout
* ``Intervals.xlsx``     background/signal intervals in window seconds, the
  same for every sample (arguments), oxide wt% columns left blank for the
  autofill / the user
* ``DriftSelections.xlsx`` / ``CalibrationSelections.xlsx`` defaults, optionally
  mirrored per element from a previous run's selection files
"""
from __future__ import annotations

import glob
import io
import os
import re
from datetime import datetime
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from . import io_xlsx
from .names import mangle_columns
from .core import INTERVAL_COLUMNS

TS_FORMATS = ('%m/%d/%Y %I:%M:%S %p', '%m/%d/%Y %H:%M:%S')


def parse_dwell_line(line: str, columns: Sequence[str]) -> Dict[str, float]:
    """The ',dwell time=0.01;xcal factor=...' line -> {isotope: dwell seconds}."""
    out = {}
    cells = line.split(',')
    for col, cell in zip(columns, cells[1:]):
        m = re.search(r'dwell time=([0-9.eE+-]+)', cell)
        if m:
            try:
                out[str(col)] = float(m.group(1))
            except ValueError:
                pass
    return out


def parse_file(path: str, want_dwell: bool = False):
    """-> (sample name, acquisition timestamp, DataFrame with Time + isotopes[, dwell dict])."""
    with open(path, encoding='utf-8', errors='replace') as f:
        lines = f.read().splitlines()
    first = lines[0].rstrip(';').strip()
    name, _, ts = first.partition(':')
    ts = ts.strip()
    stamp = None
    for fmt in TS_FORMATS:
        try:
            stamp = datetime.strptime(ts, fmt)
            break
        except ValueError:
            pass
    if stamp is None:
        raise ValueError(f'{path}: cannot parse timestamp "{ts}"')
    hdr = next(i for i, l in enumerate(lines) if l.startswith('Time,'))
    body = [lines[hdr]] + [l for l in lines[hdr + 1:] if l and not l.startswith(',dwell')]
    df = pd.read_csv(io.StringIO('\n'.join(body)), index_col=False)
    df = df.loc[:, [c for c in df.columns if not str(c).startswith('Unnamed')]]
    df = df.apply(pd.to_numeric, errors='coerce')
    if want_dwell:
        dwell_line = next((l for l in lines[hdr + 1:] if l.startswith(',dwell')), '')
        cols = [c for c in lines[hdr].split(',')[1:] if c]
        return name.strip(), stamp, df, parse_dwell_line(dwell_line, cols)
    return name.strip(), stamp, df


def file_number(path: str) -> int:
    m = re.search(r'_(\d+)\.csv$', os.path.basename(path))
    return int(m.group(1)) if m else -1


def list_files(src_dir: str) -> List[str]:
    files = glob.glob(os.path.join(src_dir, '*.csv'))
    return sorted(files, key=file_number)


def default_cal_prefs(sig_vars: Sequence[str], template: Optional[pd.DataFrame],
                      norm_default: str = 'Al') -> pd.DataFrame:
    """Per-analyte NormElement / ForceInterceptZero / StandardSet, mirrored by
    element symbol from ``template`` (a previous CalibrationSelections) when
    the symbol exists there, otherwise (norm_default, True, GeoRem).
    Exclusions are never carried over."""
    def symbol(v):
        return re.sub(r'^x?\d+', '', str(v))
    tmpl = {}
    if template is not None:
        for row, r in template.iterrows():
            tmpl[symbol(row)] = r
    rows = []
    for v in sig_vars:
        r = tmpl.get(symbol(v))
        if r is not None:
            rows.append({'NormElement': str(r['NormElement']), 'ForceInterceptZero': bool(r['ForceInterceptZero']),
                         'StandardSet': str(r['StandardSet']), 'CalExclusions': 'None',
                         'Weighting': str(r['Weighting']) if 'Weighting' in r.index else 'ivw',
                         'Calibrants': str(r['Calibrants']) if 'Calibrants' in r.index else 'BHV BCR BIR'})
        else:
            rows.append({'NormElement': norm_default, 'ForceInterceptZero': True,
                         'StandardSet': 'GeoRem', 'CalExclusions': 'None', 'Weighting': 'ivw',
                         'Calibrants': 'BHV BCR BIR'})
    return pd.DataFrame(rows, index=list(sig_vars), dtype=object)


def build_data_folder(src_dir: str, out_dir: str, name: str,
                      back1: Tuple[float, float] = (2.0, 24.0),
                      signal: Tuple[float, float] = (29.5, 45.5),
                      back2: Optional[Tuple[float, float]] = None,
                      start_thresh: float = 1000.0, start_iso: str = 'x25Mg',
                      cal_template: Optional[pd.DataFrame] = None,
                      drift_std: str = 'ALL', drift_interp: str = 'auto',
                      norm_default: str = 'Al') -> Dict[str, str]:
    os.makedirs(out_dir, exist_ok=True)
    files = list_files(src_dir)
    if not files:
        raise FileNotFoundError(f'no CSV files in {src_dir}')
    parsed4 = [parse_file(f, want_dwell=True) for f in files]
    parsed = [p[:3] for p in parsed4]
    t0 = min(p[1] for p in parsed)
    dwell_raw: Dict[str, float] = {}
    for p in parsed4:
        for k, v in p[3].items():
            dwell_raw.setdefault(k, v)

    names, frames, windows, meta = [], [], [], []
    row0 = 0
    for path, (nm, stamp, df) in zip(files, parsed):
        df = df.copy()
        df.columns = mangle_columns(df.columns)
        offset = (stamp - t0).total_seconds()
        local_t = df['Time'].to_numpy(float)
        df['Time'] = offset + local_t
        n = len(df)
        on = np.flatnonzero(df[start_iso].to_numpy(float) >= start_thresh) if start_iso in df.columns else np.array([], int)
        on_row = int(on[0]) if on.size else 0
        names.append(nm)
        frames.append(df)
        windows.append({'start_row': row0 + 1, 'stop_row': row0 + n, 'laser_on_row': row0 + on_row + 1})
        meta.append({'File': os.path.basename(path), 'FileNumber': file_number(path), 'SampleName': nm,
                     'Timestamp': stamp.strftime('%Y-%m-%d %H:%M:%S'), 'OffsetSeconds': offset,
                     'Rows': n, 'LaserOnLocalSeconds': float(local_t[on_row]) if on.size else np.nan})
        row0 += n
    batch = pd.concat(frames, ignore_index=True)

    data_path = os.path.join(out_dir, f'{name}_data.xlsx')
    io_xlsx.write_table(data_path, batch, 'Sheet1')
    ro_path = os.path.join(out_dir, 'RunOrder.xlsx')
    io_xlsx.write_run_order(ro_path, names)
    sw_path = os.path.join(out_dir, 'SampleWindows.xlsx')
    from .names import clean_sample_names
    disp = clean_sample_names(names)
    io_xlsx.write_row_table(sw_path, pd.DataFrame(windows, index=disp), 'Sheet1')
    meta_path = os.path.join(out_dir, 'SourceFiles.xlsx')
    io_xlsx.write_table(meta_path, pd.DataFrame(meta), 'Sheet1')
    dwell = {m: dwell_raw[c] for c, m in zip(dwell_raw.keys(), mangle_columns(list(dwell_raw.keys())))}
    dwell_path = os.path.join(out_dir, 'DwellTimes.xlsx')
    if dwell:
        io_xlsx.write_dwell_times(dwell_path, dwell)

    iv = pd.DataFrame(np.nan, index=disp, columns=INTERVAL_COLUMNS)
    iv['back_start1'], iv['back_stop1'] = back1
    iv['signal_start'], iv['signal_stop'] = signal
    if back2 is not None:
        iv['back_start2'], iv['back_stop2'] = back2
    iv_path = os.path.join(out_dir, 'Intervals.xlsx')
    io_xlsx.write_intervals(iv_path, iv)

    sig_vars = [c for c in batch.columns if c != 'Time']
    dp = pd.DataFrame({'Standard': drift_std, 'Interpolation': drift_interp, 'Exclusions': ''},
                      index=sig_vars, dtype=object)
    dp_path = os.path.join(out_dir, 'DriftSelections.xlsx')
    io_xlsx.write_drift_selections(dp_path, dp)
    cp = default_cal_prefs(sig_vars, cal_template, norm_default)
    cp_path = os.path.join(out_dir, 'CalibrationSelections.xlsx')
    io_xlsx.write_calibration_selections(cp_path, cp)
    return {'data': data_path, 'run_order': ro_path, 'windows': sw_path, 'intervals': iv_path,
            'drift': dp_path, 'cal': cp_path, 'source_files': meta_path, 'n_samples': len(names)}
