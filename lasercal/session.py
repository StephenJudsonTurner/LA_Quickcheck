"""ReductionSession: the state machine behind LaserCalTool.

It reproduces, in order, the MATLAB start-up sequence (load data, RunOrder,
find starts, clean names, build windows, load/create Intervals.xlsx and
autofill oxide wt%), then exposes the three button actions:

* run_drift()          -> "Calculate drifts"
* calibrate()          -> "Compute"
* export_drift_corrected() / export_calibration() -> the two "Export" buttons

plus the two "Save ... preferences" actions.  The MATLAB code falls back to
the *current UI dropdown values* whenever a preference is missing; those
values live in ``UIState`` so the GUI and the headless runner share one code
path.
"""
from __future__ import annotations

import os
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence

import numpy as np
import pandas as pd

from . import core, io_xlsx
from .names import clean_sample_names

STD_FILE_NAME = 'LAICPMS_Standard_Values.xlsx'


@dataclass
class UIState:
    """Current values of the dropdowns / checkboxes the MATLAB code reads as
    fall-backs.  Defaults equal the MATLAB initial UI state."""
    drift_std: str = 'BHV'
    drift_interp: str = 'pchip'
    cal_norm: Optional[str] = None   # None -> first available choice (showChoices{1})
    force_zero: bool = True
    harvard: bool = False


class IntervalsCreated(Exception):
    """Raised when Intervals.xlsx did not exist and was created (MATLAB shows
    an alert and exits)."""


class ReductionSession:
    def __init__(self, data_file: str, std_file: Optional[str] = None,
                 std_vals: Optional[pd.DataFrame] = None, ui: Optional[UIState] = None,
                 program_dir: Optional[str] = None):
        self.data_file = os.path.abspath(data_file)
        self.folder = os.path.dirname(self.data_file)
        self.program_dir = program_dir or os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        self.std_file = std_file or os.path.join(self.program_dir, STD_FILE_NAME)
        self._std_vals_override = std_vals
        self.ui = ui or UIState()

        # populated by load()
        self.batch: Optional[pd.DataFrame] = None
        self.sample_names: List[str] = []
        self.display_names: List[str] = []
        self.first_in: np.ndarray = np.array([], int)
        self.first_in_times: np.ndarray = np.array([])
        self.windows: List[pd.DataFrame] = []
        self.intervals: Optional[pd.DataFrame] = None
        self.sig_vars: List[str] = []
        self.norm_map = {}
        self.show_choices: List[str] = []
        self.sets: Dict[str, np.ndarray] = {}

        # drift state
        self.av_total: Optional[np.ndarray] = None
        self.av_back: Optional[np.ndarray] = None
        self.drift_grid: Optional[np.ndarray] = None
        self.lt_corr: Optional[np.ndarray] = None
        self.sem_pct: Optional[np.ndarray] = None

        # calibration state
        self.cal_elements: List[core.CalElement] = []
        self.final_result: Optional[np.ndarray] = None
        self.std_vals: Optional[pd.DataFrame] = None
        self.cal_excl_memory: Dict[str, List[str]] = {}   # CalExclByAnalyte (names only)

    # ------------------------------------------------------------------ paths
    @property
    def intervals_path(self): return os.path.join(self.folder, 'Intervals.xlsx')
    @property
    def run_order_path(self): return os.path.join(self.folder, 'RunOrder.xlsx')
    @property
    def drift_prefs_path(self): return os.path.join(self.folder, 'DriftSelections.xlsx')
    @property
    def cal_prefs_path(self): return os.path.join(self.folder, 'CalibrationSelections.xlsx')
    @property
    def sample_windows_path(self): return os.path.join(self.folder, 'SampleWindows.xlsx')

    @property
    def n_samples(self): return len(self.display_names)

    @property
    def cal_norm_default(self) -> str:
        if self.ui.cal_norm and self.ui.cal_norm in self.norm_map:
            return self.ui.cal_norm
        return self.show_choices[0] if self.show_choices else 'Ca'

    # ------------------------------------------------------------------ startup
    def load(self):
        """Everything LaserCalTool does before building the figure."""
        self.batch = io_xlsx.read_raw_data(self.data_file)
        for req in ('Time', 'x25Mg'):
            if req not in self.batch.columns:
                raise ValueError('Data file must include "Time" and "x25Mg" columns.')
        self.sample_names = io_xlsx.read_run_order(self.run_order_path)
        n = len(self.sample_names)
        self.display_names = clean_sample_names(self.sample_names)
        if os.path.isfile(self.sample_windows_path):
            # Explicit windows (per-analysis file exports): each sample is the row
            # span given in SampleWindows.xlsx; "start time" = the laser-on row.
            sw = io_xlsx.read_sample_windows(self.sample_windows_path)
            if len(sw) != n:
                raise ValueError(f'SampleWindows.xlsx has {len(sw)} rows but RunOrder.xlsx has {n}.')
            a = sw['start_row'].to_numpy(int) - 1
            b = sw['stop_row'].to_numpy(int) - 1
            self.first_in = sw['laser_on_row'].to_numpy(int) - 1
            self.windows = []
            for s, e in zip(a, b):
                T = self.batch.iloc[s:e + 1].reset_index(drop=True).copy()
                T['Time'] = T['Time'] - T['Time'].iloc[0]
                self.windows.append(T)
            self.window_mode = 'explicit'
        else:
            self.first_in = core.find_sample_starts(self.batch['x25Mg'].to_numpy(float), n)
            self.windows = core.build_windows(self.batch, self.first_in)
            self.window_mode = 'threshold'
        self.first_in_times = self.batch['Time'].to_numpy(float)[self.first_in]

        all_vars = list(self.windows[0].columns)
        self.sig_vars = [v for v in all_vars if v != 'Time']
        if not self.sig_vars:
            raise ValueError('No signal columns found in data.')
        self.norm_map = core.build_norm_map(self.sig_vars)
        self.show_choices = core.norm_choices(self.sig_vars, self.norm_map)
        self.sets = core.standard_idx_sets(self.display_names)

        if os.path.isfile(self.intervals_path):
            self.intervals = io_xlsx.read_intervals(self.intervals_path)
        else:
            iv = pd.DataFrame(np.full((n, len(core.INTERVAL_COLUMNS)), np.nan),
                              columns=core.INTERVAL_COLUMNS, index=self.display_names)
            io_xlsx.write_intervals(self.intervals_path, iv)
            raise IntervalsCreated('Intervals.xlsx created in selected folder. Please fill it out and re-run.')

        self.intervals = core.autofill_oxides(self.intervals, self.sets)
        io_xlsx.write_intervals(self.intervals_path, self.intervals)
        return self

    # ------------------------------------------------------------------ preferences
    def load_or_init_drift_prefs(self) -> pd.DataFrame:
        """loadOrInitDriftPrefs: rows = sig_vars, cols Standard/Interpolation/Exclusions."""
        cols = ['Standard', 'Interpolation', 'Exclusions']
        fp = self.drift_prefs_path
        write_back = False
        prefs = None
        if os.path.isfile(fp):
            try:
                prefs = io_xlsx.read_row_table(fp)
                prefs = prefs[[c for c in cols if c in prefs.columns]]
            except Exception:
                prefs = None
        if prefs is None:
            prefs = pd.DataFrame({'Standard': 'BHV', 'Interpolation': 'pchip', 'Exclusions': ''},
                                 index=self.sig_vars, dtype=object)
            write_back = not os.path.isfile(fp)
        for c in cols:
            if c not in prefs.columns:
                prefs[c] = ''
        missing = [v for v in self.sig_vars if v not in prefs.index]
        if missing:
            add = pd.DataFrame({'Standard': 'BHV', 'Interpolation': 'pchip', 'Exclusions': ''},
                               index=missing, dtype=object)
            prefs = pd.concat([prefs, add])
            write_back = True
        prefs = prefs.loc[self.sig_vars, cols].astype(object)
        if write_back:
            try:
                io_xlsx.write_drift_selections(fp, prefs)
            except Exception:
                pass
        return prefs

    def load_or_init_cal_prefs(self) -> pd.DataFrame:
        """loadOrInitCalibrationPrefs (never writes the file)."""
        cols = ['NormElement', 'ForceInterceptZero', 'StandardSet', 'CalExclusions']
        fp = self.cal_prefs_path
        default_key = self.cal_norm_default
        prefs = None
        if os.path.isfile(fp):
            try:
                prefs = io_xlsx.read_row_table(fp)
            except Exception:
                prefs = None
        if prefs is None:
            prefs = pd.DataFrame(index=self.sig_vars, dtype=object)
        n = len(prefs)
        if 'NormElement' not in prefs.columns:
            prefs['NormElement'] = [default_key] * n
        if 'ForceInterceptZero' not in prefs.columns:
            prefs['ForceInterceptZero'] = [True] * n
        if 'StandardSet' not in prefs.columns:
            prefs['StandardSet'] = ['GeoRem'] * n
        if 'CalExclusions' not in prefs.columns:
            prefs['CalExclusions'] = ['None'] * n
        missing = [v for v in self.sig_vars if v not in prefs.index]
        if missing:
            add = pd.DataFrame({'NormElement': default_key, 'ForceInterceptZero': True,
                                'StandardSet': 'GeoRem', 'CalExclusions': 'None'}, index=missing, dtype=object)
            prefs = pd.concat([prefs, add])
        return prefs.loc[self.sig_vars, cols].astype(object)

    def save_drift_pref(self, el: str, std: str, interp: str, excl_idx: Sequence[int]) -> bool:
        """saveElementPreference. excl_idx are 1-based sample numbers currently
        marked in the plot. Returns True if the standard changed (exclusions cleared)."""
        prefs = self.load_or_init_drift_prefs()
        if el not in prefs.index:
            prefs.loc[el] = [std, interp, '']
        prev = str(prefs.loc[el, 'Standard']) if not core.is_blank(prefs.loc[el, 'Standard']) else ''
        std_changed = prev.lower() != std.lower()
        if std_changed:
            exstr = ''
        else:
            ex = sorted(set(int(i) for i in excl_idx))
            exstr = ' '.join(str(i) for i in ex) if ex else ''
        prefs.loc[el, 'Standard'] = std
        prefs.loc[el, 'Interpolation'] = interp
        prefs.loc[el, 'Exclusions'] = exstr
        io_xlsx.write_drift_selections(self.drift_prefs_path, prefs)
        return std_changed

    def save_cal_pref(self, el: str, norm_key: str, force_zero: bool, harvard: bool) -> str:
        """saveCalibrationPreference. Uses the in-memory exclusions for `el`."""
        prefs = self.load_or_init_cal_prefs()
        names = [n for n in self.cal_excl_memory.get(el, []) if n and n != 'None']
        exstr = ' '.join(names) if names else 'None'
        prefs.loc[el, 'NormElement'] = norm_key
        prefs.loc[el, 'ForceInterceptZero'] = bool(force_zero)
        prefs.loc[el, 'StandardSet'] = 'Harvard' if harvard else 'GeoRem'
        prefs.loc[el, 'CalExclusions'] = exstr
        io_xlsx.write_calibration_selections(self.cal_prefs_path, prefs)
        return exstr

    # ------------------------------------------------------------------ drift
    def run_drift(self):
        """runDrift(): background-corrected averages, drift grid, corrected data
        and the SEM% matrix used for error bars."""
        prefs = self.load_or_init_drift_prefs()
        self.av_total, self.av_back = core.compute_averaged_signals(self.windows, self.intervals, self.sig_vars)
        try:
            cal_prefs = self.load_or_init_cal_prefs()
            self.sem_pct = core.compute_sem_backcorr_normalized(self.windows, self.intervals, self.sig_vars,
                                                                 self.norm_map, cal_prefs)
        except Exception:
            self.sem_pct = None
        self.drift_grid = core.compute_drift_grid(self.av_total, self.first_in_times, self.sig_vars, prefs,
                                                  self.sets, self.ui.drift_std, self.ui.drift_interp)
        self.lt_corr = self.av_total * (1.0 / self.drift_grid)
        return self.lt_corr

    def drift_diag(self, el: str) -> core.DriftElementDiag:
        prefs = self.load_or_init_drift_prefs()
        return core.drift_diagnostics(self.av_total, self.first_in_times, self.sig_vars, prefs, self.sets, el,
                                      self.ui.drift_std, self.ui.drift_interp)

    def export_drift_corrected(self, path: Optional[str] = None) -> str:
        if self.lt_corr is None:
            raise RuntimeError('Run drift correction first.')
        path = path or os.path.join(self.folder, 'DriftCorrected_All.xlsx')
        io_xlsx.write_drift_corrected(path, self.lt_corr, self.sig_vars, self.display_names,
                                      self.norm_map['Ca']['iso'])
        return path

    # ------------------------------------------------------------------ calibration
    def load_std_values(self) -> pd.DataFrame:
        if self._std_vals_override is not None:
            return self._std_vals_override
        if not os.path.isfile(self.std_file):
            raise FileNotFoundError(f'Standards file not found: {self.std_file}')
        return io_xlsx.read_standards(self.std_file)

    def calibrate(self):
        """calibrateAll(): fits and FinalResult."""
        if self.lt_corr is None:
            raise RuntimeError('Please run Drift Correction first.')
        self.std_vals = self.load_std_values()
        cal_prefs = self.load_or_init_cal_prefs()
        self.cal_elements, self.final_result = core.calibrate_all(
            self.lt_corr, self.sig_vars, self.display_names, self.sets, self.std_vals, cal_prefs,
            self.norm_map, self.intervals, self.cal_norm_default, self.ui.force_zero, self.ui.harvard,
            self.cal_excl_memory)
        return self.final_result

    def summary_rows(self, analyte: str):
        if self.final_result is None:
            return []
        set_name = core.cal_set_name(self.load_or_init_cal_prefs(), analyte, self.ui.harvard)
        return core.secondary_summary(self.display_names, analyte, self.sig_vars, self.final_result,
                                      self.std_vals, set_name)

    def compute_sem_for_export(self) -> Optional[np.ndarray]:
        try:
            cal_prefs = self.load_or_init_cal_prefs()
            return core.compute_sem_backcorr_normalized(self.windows, self.intervals, self.sig_vars,
                                                        self.norm_map, cal_prefs)
        except Exception:
            return None

    def export_calibration(self, path: Optional[str] = None) -> str:
        if self.final_result is None:
            raise RuntimeError('Run calibration first.')
        path = path or os.path.join(self.folder, 'ReducedDataExport.xlsx')
        sem = self.compute_sem_for_export()
        io_xlsx.write_reduced_export(path, self.final_result, self.sig_vars, self.display_names,
                                     self.first_in_times, sem, self.std_vals, self.intervals,
                                     self.load_or_init_drift_prefs(), self.load_or_init_cal_prefs())
        return path
