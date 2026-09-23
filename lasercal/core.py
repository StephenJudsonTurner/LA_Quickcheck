"""Numerical core of the reduction, mirroring LaserCalTool_Beta_17.m.

Every function here is a straight transcription of one block or nested
function of the MATLAB file (named in each docstring).  Indices are 0-based
inside Python; the only place 1-based sample numbers appear is in the
Exclusions strings of DriftSelections.xlsx, which are converted at the
boundary (see ``parse_exclusions``).
"""
from __future__ import annotations

import re
import warnings
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from .interp import interp1_extrap, polyfit_polyval

EPS = np.finfo(float).eps  # MATLAB eps

# ----------------------------------------------------------------------------
# Tunables (MATLAB "Tunables" block)
# ----------------------------------------------------------------------------
GAP = 80             # rows between sample starts
SUB = 50             # rows before/after [start : start+Gap]
START_THRESH = 1000  # single threshold on x25Mg for all samples

INTERVAL_COLUMNS = ['back_start1', 'back_stop1', 'signal_start', 'signal_stop',
                    'back_start2', 'back_stop2', 'CaO (wt%)', 'SiO2 (wt%)', 'Al2O3 (wt%)']

# Autofill defaults for known standards (MATLAB: fillDefault block).
# Tuple order is (CaO, Al2O3, SiO2) exactly as the MATLAB code assigns them.
OXIDE_DEFAULTS = {
    'BHV':     (11.14, 13.6, 49.3),
    'BIR':     (13.24, 15.5, 47.5),
    'BCR':     (6.95, 13.4, 54.4),
    'VE32':    (10.86, 13.6, 50.82),
    'GSD':     (7.2, 13.4, 53.2),
    'GSE':     (7.2, 13.4, 53.2),  # GSE-1G shares the GSD-1G basalt matrix (see standards table)
    'GOR_128': (6.24, 10, 40),     # placeholders in the MATLAB code
    'STH':     (5.28, 17.8, 63.7),
}

# Standards recognised by name (lower-cased substring of the sample name) and the
# row tag used to find them in LAICPMS_Standard_Values.xlsx.
STANDARD_TAGS = {'VE32': 've32', 'BHV': 'bhv', 'BCR': 'bcr', 'BIR': 'bir',
                 'GSD': 'gsd', 'GSE': 'gse', 'STH': 'sth', 'GOR_128': 'gor_128'}
STANDARD_ROW_KEY = {'VE32': 'VE32', 'GSD': 'GSD1G', 'GSE': 'GSE1G', 'STH': 'STHS', 'GOR_128': 'GOR-128'}
PRIMARY_STDS = ('BHV', 'BCR', 'BIR')

# Drift: a per-element "Standard" of ALL pools every recognised standard
# (count-weighted); a single key gives the legacy single-standard model.
VALID_DRIFT_STDS = ('ALL', 'BHV', 'BCR', 'GSD', 'BIR', 'GSE', 'STH')
VALID_METHODS = ('none', 'auto', 'linear', 'pchip', 'makima', 'spline', 'poly2', 'poly3')

# Counting statistics
DWELL_S_DEFAULT = 0.01      # dwell per isotope per sweep (s) when no DwellTimes.xlsx is present
ERR_FLOOR_PCT = 0.5         # added in quadrature to the counting error for weights
DRIFT_SE_SWITCH_PCT = 1.5   # 'auto' drift: interpolate when the bracket-mean error is below this, else poly2
MAX_CAL_ERR_PCT = 10.0      # calibrant points with a larger counting error are dropped from the fit

# Calibration weighting schemes (CalibrationSelections 'Weighting')
VALID_WEIGHTINGS = ('ivw', 'logmean', 'ols0', 'ols')

NORM_OXIDE_COL = {'Ca': 'CaO (wt%)', 'Al': 'Al2O3 (wt%)', 'Si': 'SiO2 (wt%)'}
NORM_FALLBACK_ISO = {'Ca': 'x43Ca', 'Al': 'x27Al', 'Si': 'x29Si'}


# ----------------------------------------------------------------------------
# Small MATLAB-semantics helpers
# ----------------------------------------------------------------------------
def nanmean_rows(a: np.ndarray) -> np.ndarray:
    """mean(A,1,'omitnan') for a 2-D array; empty input -> NaN row."""
    a = np.asarray(a, float)
    if a.shape[0] == 0:
        return np.full(a.shape[1], np.nan)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        return np.nanmean(a, axis=0)


def nanmean(v) -> float:
    v = np.asarray(v, float).ravel()
    if v.size == 0:
        return np.nan
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        return float(np.nanmean(v))


def nanstd(v) -> float:
    """std(v,'omitnan') - sample standard deviation (N-1)."""
    v = np.asarray(v, float).ravel()
    v = v[np.isfinite(v)]
    if v.size == 0:
        return np.nan
    if v.size == 1:
        return 0.0
    return float(np.std(v, ddof=1))


def is_blank(v) -> bool:
    """True for the values MATLAB would treat as an empty string / missing."""
    if v is None:
        return True
    if isinstance(v, float) and np.isnan(v):
        return True
    s = str(v).strip()
    return s == '' or s.lower() in ('nan', '<missing>')


def parse_exclusions(s) -> np.ndarray:
    """parseExclusionsString: integers found in a space/comma list (1-based).

    Non-string cells are converted the way ``string(...)`` would in MATLAB;
    the regex ``[-]?\\d+`` then pulls out integers, which are unique-sorted.
    """
    if s is None:
        return np.array([], int)
    if isinstance(s, (float, np.floating)):
        if np.isnan(s):
            return np.array([], int)
        s = str(int(s)) if float(s).is_integer() else repr(float(s))
    s = str(s)
    toks = re.findall(r"-?\d+", s)
    if not toks:
        return np.array([], int)
    return np.unique(np.array([int(t) for t in toks]))


def detect_isotope_for_element(var_list: Sequence[str], elem: str, fallback: str) -> str:
    """detectIsotopeForElement: first variable ending with the element symbol."""
    for v in var_list:
        if v.endswith(elem):
            return v
    return fallback


def build_norm_map(sig_vars: Sequence[str]) -> Dict[str, Dict[str, str]]:
    return {k: {'iso': detect_isotope_for_element(sig_vars, k, NORM_FALLBACK_ISO[k]),
                'oxideCol': NORM_OXIDE_COL[k]} for k in ('Ca', 'Al', 'Si')}


def norm_choices(sig_vars: Sequence[str], norm_map) -> List[str]:
    """showChoices: only the normalizers whose isotope exists in the data."""
    return [k for k in ('Ca', 'Al', 'Si') if norm_map[k]['iso'] in sig_vars]


# ----------------------------------------------------------------------------
# Sample starts and windows
# ----------------------------------------------------------------------------
def find_sample_starts(x25mg: np.ndarray, n_samples: int,
                       gap: int = GAP, thresh: float = START_THRESH) -> np.ndarray:
    """'Locate per-sample starts' block. Returns 0-based row indices."""
    mg = np.asarray(x25mg, float)
    n = mg.size
    hits = np.flatnonzero(mg >= thresh)
    if hits.size == 0:
        raise ValueError(f"Could not find any start >= {thresh:g}.")
    starts = [int(hits[0])]
    for i in range(2, n_samples + 1):
        start_count = starts[-1] + gap            # MATLAB: FirstIn(end)+Gap
        start_count = min(start_count, n - 1)     # MATLAB: min(StartCount, height)
        sub = mg[start_count:]
        k = np.flatnonzero(sub >= thresh)
        if k.size == 0:
            raise ValueError(f"Could not find subsequent start >= {thresh:g} for sample #{i}.")
        starts.append(start_count + int(k[0]))
    return np.array(starts, int)


def build_windows(batch: pd.DataFrame, starts: np.ndarray,
                  gap: int = GAP, sub: int = SUB) -> List[pd.DataFrame]:
    """'Build MyData windows (raw)'. Each window's Time is re-zeroed."""
    n = len(batch)
    wins = []
    for s in starts:
        a = max(s - sub, 0)
        b = min(s + gap + sub, n - 1)
        T = batch.iloc[a:b + 1].reset_index(drop=True).copy()
        if 'Time' in T.columns and len(T):
            T['Time'] = T['Time'] - T['Time'].iloc[0]
        wins.append(T)
    return wins


# ----------------------------------------------------------------------------
# Standards
# ----------------------------------------------------------------------------
def standard_idx_sets(display_names: Sequence[str]) -> Dict[str, np.ndarray]:
    """getStandardIdxSets: 0-based indices of samples whose lower-cased name
    contains the tag."""
    names = [str(n).lower() for n in display_names]
    return {k: np.array([i for i, n in enumerate(names) if t in n], int) for k, t in STANDARD_TAGS.items()}


def all_standard_idx(sets: Dict[str, np.ndarray]) -> np.ndarray:
    """Sorted union of every recognised standard's indices."""
    return np.unique(np.concatenate([np.asarray(v, int) for v in sets.values()] + [np.array([], int)]))


def standard_brackets(sets: Dict[str, np.ndarray]) -> List[np.ndarray]:
    """Brackets = runs of consecutive standard analyses, split where a standard
    that is already in the current run appears again (two brackets run back to
    back stay two brackets)."""
    union = all_standard_idx(sets)
    key_of = {int(i): k for k, v in sets.items() for i in v}
    groups: List[List[int]] = []
    cur: List[int] = []
    seen: set = set()
    for i in union:
        i = int(i)
        k = key_of.get(i)
        if cur and (i != cur[-1] + 1 or k in seen):
            groups.append(cur); cur = []; seen = set()
        cur.append(i); seen.add(k)
    if cur:
        groups.append(cur)
    return [np.array(g, int) for g in groups]


def standard_key_of(i: int, sets: Dict[str, np.ndarray]) -> Optional[str]:
    for k, v in sets.items():
        if i in v:
            return k
    return None


def autofill_oxides(intervals: pd.DataFrame, sets: Dict[str, np.ndarray]) -> pd.DataFrame:
    """Fill NaN CaO/Al2O3/SiO2 for known standards (positional rows)."""
    iv = intervals.copy()
    for key, (cao, al2o3, sio2) in OXIDE_DEFAULTS.items():
        idx = sets.get(key, np.array([], int))
        if idx.size == 0:
            continue
        for col, val in (('CaO (wt%)', cao), ('Al2O3 (wt%)', al2o3), ('SiO2 (wt%)', sio2)):
            if col not in iv.columns:
                continue
            colvals = pd.to_numeric(iv[col], errors='coerce').to_numpy(dtype=float, copy=True)
            m = np.isnan(colvals[idx])
            colvals[idx[m]] = val
            iv[col] = colvals
    return iv


# ----------------------------------------------------------------------------
# Background-corrected averages
# ----------------------------------------------------------------------------
def _iv(intervals: pd.DataFrame, row: int, col: str) -> float:
    if col not in intervals.columns or row >= len(intervals):
        return np.nan
    v = intervals[col].iloc[row]
    try:
        return float(v)
    except (TypeError, ValueError):
        return np.nan


@dataclass
class SignalStats:
    """Window means plus what counting statistics need."""
    av_total: np.ndarray      # background-corrected mean cps (n x m)
    av_back: np.ndarray       # background mean cps
    av_gross: np.ndarray      # signal-window mean cps (not background-corrected)
    n_signal: np.ndarray      # sweeps in the signal window (n,)
    n_back: np.ndarray        # sweeps in the background window(s) (n,)


def compute_signal_stats(windows: List[pd.DataFrame], intervals: pd.DataFrame,
                         sig_order: Sequence[str]) -> SignalStats:
    """computeAveragedSignals_raw plus gross means and sweep counts."""
    av_total, av_back, av_gross, n_sig, n_back = _averaged_signals(windows, intervals, sig_order)
    return SignalStats(av_total, av_back, av_gross, n_sig, n_back)


def counting_error_pct(stats: SignalStats, dwell_s: np.ndarray) -> np.ndarray:
    """Relative Poisson error (%) of each background-corrected mean.

    var(net cps) = gross/(n_signal*dwell) + back/(n_back*dwell); dwell_s is one
    value per isotope (s).  Non-positive nets give +inf (no weight).
    """
    dwell = np.asarray(dwell_s, float)[None, :]
    with np.errstate(divide='ignore', invalid='ignore'):
        var = (np.fmax(stats.av_gross, 0) / (stats.n_signal[:, None] * dwell)
               + np.fmax(stats.av_back, 0) / (np.fmax(stats.n_back, 1)[:, None] * dwell))
        err = 100 * np.sqrt(var) / stats.av_total
    err = np.where(np.isfinite(err) & (stats.av_total > 0), err, np.inf)
    return err


def ratio_error_pct(cnt_err: np.ndarray, j: int, ni: int) -> np.ndarray:
    """Counting error (%) of the ratio isotope j / isotope ni."""
    return np.sqrt(cnt_err[:, j] ** 2 + cnt_err[:, ni] ** 2)


def compute_averaged_signals(windows: List[pd.DataFrame], intervals: pd.DataFrame,
                             sig_order: Sequence[str]) -> Tuple[np.ndarray, np.ndarray]:
    """computeAveragedSignals_raw -> (AvTotalCorrected, AvBack)."""
    av_total, av_back, _, _, _ = _averaged_signals(windows, intervals, sig_order)
    return av_total, av_back


def _averaged_signals(windows, intervals, sig_order):
    n = len(windows)
    m = len(sig_order)
    av_back = np.full((n, m), np.nan)
    av_int1 = np.full((n, m), np.nan)
    av_int2 = np.full((n, m), np.nan)
    n_sig = np.zeros(n, int)
    n_back = np.zeros(n, int)
    for ii, T in enumerate(windows):
        t = T['Time'].to_numpy(float)
        S = T[list(sig_order)].to_numpy(float)
        s0, s1 = _iv(intervals, ii, 'signal_start'), _iv(intervals, ii, 'signal_stop')
        mid = (s0 + s1) / 2
        b1, e1 = _iv(intervals, ii, 'back_start1'), _iv(intervals, ii, 'back_stop1')
        b2, e2 = _iv(intervals, ii, 'back_start2'), _iv(intervals, ii, 'back_stop2')
        with np.errstate(invalid='ignore'):
            m_b1 = (~np.isnan(b1)) & (~np.isnan(e1)) & (t >= b1) & (t < e1)
            m_b2 = (~np.isnan(b2)) & (~np.isnan(e2)) & (t >= b2) & (t < e2)
            m_s1 = (t >= s0) & (t < mid)
            m_s2 = (t >= mid) & (t < s1)
        B = S[m_b1]
        if m_b2.any():
            B = np.vstack([B, S[m_b2]])
        av_back[ii] = nanmean_rows(B)
        av_int1[ii] = nanmean_rows(S[m_s1])
        av_int2[ii] = nanmean_rows(S[m_s2])
        n_sig[ii] = int(m_s1.sum() + m_s2.sum())
        n_back[ii] = int(B.shape[0])
    av_gross = (av_int1 + av_int2) / 2
    return av_gross - av_back, av_back, av_gross, n_sig, n_back


# ----------------------------------------------------------------------------
# Drift
# ----------------------------------------------------------------------------
def group_consecutive(idxs: np.ndarray) -> List[np.ndarray]:
    """groupStandardsByCount(idxs, Inf): runs of adjacent indices."""
    idxs = np.asarray(idxs, int).ravel()
    if idxs.size == 0:
        return []
    groups, cur = [], [int(idxs[0])]
    for k in range(1, idxs.size):
        if idxs[k] == idxs[k - 1] + 1:
            cur.append(int(idxs[k]))
        else:
            groups.append(np.array(cur, int))
            cur = [int(idxs[k])]
    groups.append(np.array(cur, int))
    return groups


@dataclass
class DriftPref:
    standard: str
    method: str
    exclusions: np.ndarray  # 1-based sample numbers


@dataclass
class DriftElementDiag:
    """What the drift tab plots for one element."""
    std_name: str
    method: str
    drift_idx_sel: np.ndarray = field(default_factory=lambda: np.array([], int))  # 0-based
    x_pts: np.ndarray = field(default_factory=lambda: np.array([]))
    y_pts: np.ndarray = field(default_factory=lambda: np.array([]))
    group_times: np.ndarray = field(default_factory=lambda: np.array([]))
    group_scaled: np.ndarray = field(default_factory=lambda: np.array([]))
    group_err_pct: np.ndarray = field(default_factory=lambda: np.array([]))
    method_used: str = ''
    ok: bool = False


def resolve_drift_pref(prefs: pd.DataFrame, el: str, ui_std: str, ui_method: str) -> DriftPref:
    """Per-element preference with the MATLAB fall-backs to the UI dropdowns."""
    std = '' if el not in prefs.index or is_blank(prefs.loc[el, 'Standard']) else str(prefs.loc[el, 'Standard']).upper()
    if std in ('POOLED', 'ALL6', '*'):
        std = 'ALL'
    if std not in VALID_DRIFT_STDS:
        std = ui_std
    meth = '' if el not in prefs.index or is_blank(prefs.loc[el, 'Interpolation']) else str(prefs.loc[el, 'Interpolation']).lower()
    if meth not in VALID_METHODS:
        meth = ui_method
    excl = parse_exclusions(prefs.loc[el, 'Exclusions']) if el in prefs.index else np.array([], int)
    return DriftPref(std, meth, excl)


def drift_curve(method: str, xg: np.ndarray, yg: np.ndarray, t_all: np.ndarray) -> np.ndarray:
    """The interpolation / fit switch of runDrift, including fall-backs."""
    method = method.lower()
    if method in ('linear', 'pchip', 'makima', 'spline'):
        try:
            return interp1_extrap(xg, yg, t_all, method)
        except Exception:
            return interp1_extrap(xg, yg, t_all, 'linear')
    if method in ('poly2', 'poly3'):
        deg = 2 if method == 'poly2' else 3
        if xg.size >= deg + 1 and np.all(np.isfinite(xg)) and np.all(np.isfinite(yg)):
            yhat = polyfit_polyval(xg, yg, t_all, deg)
            return np.where(~np.isfinite(yhat) | (yhat <= 0), EPS, yhat)
        return interp1_extrap(xg, yg, t_all, 'linear')
    return interp1_extrap(xg, yg, t_all, 'linear')


def _auto_method(n_groups: int, median_se_pct: float) -> str:
    """'auto': interpolate through the bracket means when they are precise,
    otherwise a low-order fit that averages the noise out."""
    if n_groups < 2:
        return 'none'
    if n_groups == 2:
        return 'linear'
    if median_se_pct <= DRIFT_SE_SWITCH_PCT:
        return 'pchip'
    return 'poly2' if n_groups >= 4 else 'linear'


def _weights(err_pct: np.ndarray) -> np.ndarray:
    e = np.asarray(err_pct, float)
    w = 1.0 / (e ** 2 + ERR_FLOOR_PCT ** 2)
    return np.where(np.isfinite(w), w, 0.0)


@dataclass
class PooledDrift:
    """Bracket means of one element's normalised ratio over all standards."""
    idx: np.ndarray             # 0-based standard indices used
    rel: np.ndarray             # ratio / own-standard weighted mean, per idx
    w: np.ndarray               # weights per idx
    group_times: np.ndarray
    group_vals: np.ndarray      # weighted bracket means, normalised to mean 1
    group_se_pct: np.ndarray    # 1/sqrt(sum w) per bracket (%)


def pooled_bracket_means(ratio: np.ndarray, err_pct: np.ndarray, t_all: np.ndarray,
                         sets: Dict[str, np.ndarray], std_keys: Sequence[str],
                         exclusions_1based: np.ndarray) -> Optional[PooledDrift]:
    """Count-weighted pooling of every standard in ``std_keys``.

    Each standard's ratios are divided by that standard's own weighted run
    mean (so standards of any concentration can be pooled), then the
    relative values are averaged per bracket (run of consecutive standard
    analyses) with weights 1/(err^2 + floor^2).
    """
    keys = [k for k in std_keys if sets.get(k, np.array([], int)).size > 0]
    if not keys:
        return None
    ex = np.asarray(exclusions_1based, int) - 1
    idx_all, rel_all, w_all = [], [], []
    for k in keys:
        idx = np.setdiff1d(sets[k], ex)
        r = ratio[idx]
        w = _weights(err_pct[idx])
        good = np.isfinite(r) & (r > 0) & (w > 0)
        idx, r, w = idx[good], r[good], w[good]
        if idx.size == 0 or w.sum() <= 0:
            continue
        mu = np.sum(w * r) / np.sum(w)
        if not np.isfinite(mu) or mu <= 0:
            continue
        idx_all.append(idx); rel_all.append(r / mu); w_all.append(w)
    if not idx_all:
        return None
    idx = np.concatenate(idx_all); rel = np.concatenate(rel_all); w = np.concatenate(w_all)
    order = np.argsort(idx)
    idx, rel, w = idx[order], rel[order], w[order]
    # brackets: runs of consecutive analyses among ALL recognised standards (not only the
    # ones that survived the cut), so a bracket stays one group when one member is dropped
    groups = standard_brackets(sets)
    gt, gv, gse = [], [], []
    for g in groups:
        m = np.isin(idx, g)
        if not m.any():
            continue
        sw = np.sum(w[m])
        gt.append(np.mean(t_all[g]))
        gv.append(np.sum(w[m] * rel[m]) / sw)
        gse.append(1.0 / np.sqrt(sw))          # weights are 1/%^2, so this is in %
    gt, gv, gse = map(np.asarray, (gt, gv, gse))
    den = np.mean(gv)
    if not np.isfinite(den) or den <= 0:
        return None
    return PooledDrift(idx, rel, w, gt, gv / den, gse)


def _single_standard_groups(av: np.ndarray, t_all: np.ndarray, idx: np.ndarray):
    groups = group_consecutive(idx)
    gt = np.array([np.mean(t_all[g]) for g in groups])
    gv = np.array([nanmean(av[g]) for g in groups])
    den = nanmean(gv)
    if not np.isfinite(den) or den == 0:
        return None, None
    return gt, gv / den


def compute_drift_grid(av_total: np.ndarray, first_in_times: np.ndarray, sig_order: Sequence[str],
                       prefs: pd.DataFrame, sets: Dict[str, np.ndarray],
                       ui_std: str = 'ALL', ui_method: str = 'auto',
                       cnt_err: Optional[np.ndarray] = None,
                       norm_idx: Optional[Sequence[int]] = None) -> np.ndarray:
    """Multiplicative drift factor per (analysis, element).

    * Standard = ALL (default): the element's ratio to its normaliser
      (``norm_idx[s]``) is pooled over every recognised standard, count-
      weighted (``cnt_err`` in %), and the curve goes through the bracket
      means.  The factor then applies to the *ratio* (see ``calibrate_all``).
    * Standard = one key: the legacy single-standard model on the raw signal.
    Returns a grid of 1 where no correction is possible.
    """
    n_samples, n_sig = av_total.shape
    t_all = np.asarray(first_in_times, float).ravel()
    grid = np.ones((n_samples, n_sig))
    for s, el in enumerate(sig_order):
        p = resolve_drift_pref(prefs, el, ui_std, ui_method)
        if p.method == 'none':
            continue
        if p.standard == 'ALL':
            if norm_idx is None or cnt_err is None:
                raise ValueError('pooled drift needs norm_idx and cnt_err')
            ni = int(norm_idx[s])
            if ni == s:
                continue                      # the normaliser itself carries no ratio drift
            with np.errstate(divide='ignore', invalid='ignore'):
                ratio = av_total[:, s] / av_total[:, ni]
            err = ratio_error_pct(cnt_err, s, ni)
            pooled = pooled_bracket_means(ratio, err, t_all, sets, list(STANDARD_TAGS), p.exclusions)
            if pooled is None or pooled.group_times.size < 2:
                continue
            meth = _auto_method(pooled.group_times.size, float(np.median(pooled.group_se_pct))) if p.method == 'auto' else p.method
            if meth == 'none':
                continue
            grid[:, s] = drift_curve(meth, pooled.group_times, pooled.group_vals, t_all)
        else:
            idx = sets.get(p.standard, np.array([], int))
            if p.exclusions.size:
                idx = np.setdiff1d(idx, p.exclusions - 1)   # 1-based -> 0-based
            if idx.size == 0:
                continue
            gt, gv = _single_standard_groups(av_total[:, s], t_all, idx)
            if gt is None:
                continue
            meth = _auto_method(gt.size, 0.0) if p.method == 'auto' else p.method
            if meth == 'none':
                continue
            grid[:, s] = drift_curve(meth, gt, gv, t_all)
    return grid


def drift_diagnostics(av_total: np.ndarray, first_in_times: np.ndarray, sig_order: Sequence[str],
                      prefs: pd.DataFrame, sets: Dict[str, np.ndarray], el: str,
                      ui_std: str = 'ALL', ui_method: str = 'auto',
                      cnt_err: Optional[np.ndarray] = None,
                      norm_idx: Optional[Sequence[int]] = None) -> DriftElementDiag:
    """Points shown on the 'Scaled raw drift vs fit' axes for one element."""
    t_all = np.asarray(first_in_times, float).ravel()
    p = resolve_drift_pref(prefs, el, ui_std, ui_method)
    d = DriftElementDiag(p.standard, p.method)
    sig_order = list(sig_order)
    v = sig_order.index(el) if el in sig_order else 0
    if p.standard == 'ALL':
        if norm_idx is None or cnt_err is None:
            return d
        ni = int(norm_idx[v])
        if ni == v:
            return d
        with np.errstate(divide='ignore', invalid='ignore'):
            ratio = av_total[:, v] / av_total[:, ni]
        pooled = pooled_bracket_means(ratio, ratio_error_pct(cnt_err, v, ni), t_all, sets,
                                      list(STANDARD_TAGS), p.exclusions)
        if pooled is None:
            return d
        d.drift_idx_sel = pooled.idx
        d.x_pts = t_all[pooled.idx]
        d.y_pts = pooled.rel
        d.group_times = pooled.group_times
        d.group_scaled = pooled.group_vals
        d.group_err_pct = pooled.group_se_pct
        d.method_used = _auto_method(pooled.group_times.size, float(np.median(pooled.group_se_pct))) if p.method == 'auto' else p.method
        d.ok = True
        return d
    sel = sets.get(p.standard, np.array([], int))
    if sel.size == 0:
        return d
    incl = np.setdiff1d(sel, p.exclusions - 1) if p.exclusions.size else sel
    if incl.size == 0:
        incl = sel
    gt, gv = _single_standard_groups(av_total[:, v], t_all, incl)
    if gt is None:
        return d
    den = nanmean([nanmean(av_total[g, v]) for g in group_consecutive(incl)])
    d.drift_idx_sel = sel
    d.x_pts = t_all[sel]
    d.y_pts = av_total[sel, v] / den
    d.group_times = gt
    d.group_scaled = gv
    d.method_used = _auto_method(gt.size, 0.0) if p.method == 'auto' else p.method
    d.ok = True
    return d


# ----------------------------------------------------------------------------
# SEM% of element-normalised, background-corrected signal
# ----------------------------------------------------------------------------
def sem_norm_iso_for(sig_order: Sequence[str], norm_map, cal_prefs: pd.DataFrame) -> List[str]:
    out = []
    for analyte in sig_order:
        key = 'Ca'
        if analyte in cal_prefs.index and not is_blank(cal_prefs.loc[analyte, 'NormElement']):
            key = str(cal_prefs.loc[analyte, 'NormElement'])
        if key not in norm_map:
            key = 'Ca'
        out.append(norm_map[key]['iso'])
    return out


def compute_sem_backcorr_normalized(windows: List[pd.DataFrame], intervals: pd.DataFrame,
                                    sig_order: Sequence[str], norm_map, cal_prefs: pd.DataFrame,
                                    trim_frac: float = 0.25, min_den_rel: float = 0.05,
                                    stat: str = 'sem') -> np.ndarray:
    """computeSEM_BackCorr_Normalized (smoothK = 1, i.e. no smoothing).

    stat='sem' (MATLAB): 100 * (sd/sqrt(n)) / mean of the frame ratios.
    stat='rsd'         : 100 * sd / mean  (same frames, no 1/sqrt(n)); not in MATLAB,
                         offered for plotting.
    """
    sig_order = list(sig_order)
    nS, nV = len(windows), len(sig_order)
    sem = np.full((nS, nV), np.nan)
    norm_iso_for = sem_norm_iso_for(sig_order, norm_map, cal_prefs)

    for ii, T in enumerate(windows):
        t = T['Time'].to_numpy(float)
        S = T[sig_order].to_numpy(float)
        s0, s1 = _iv(intervals, ii, 'signal_start'), _iv(intervals, ii, 'signal_stop')
        b1, e1 = _iv(intervals, ii, 'back_start1'), _iv(intervals, ii, 'back_stop1')
        b2, e2 = _iv(intervals, ii, 'back_start2'), _iv(intervals, ii, 'back_stop2')
        with np.errstate(invalid='ignore'):
            m_b1 = (~np.isnan(b1)) & (~np.isnan(e1)) & (t >= b1) & (t < e1)
            m_b2 = (~np.isnan(b2)) & (~np.isnan(e2)) & (t >= b2) & (t < e2)
            m_sig = (t >= s0) & (t < s1)
        B = S[m_b1]
        if m_b2.any():
            B = np.vstack([B, S[m_b2]])
        bmean = nanmean_rows(B)
        sig_idx = np.flatnonzero(m_sig)
        if sig_idx.size < 3:
            continue
        n_drop = int(np.floor(trim_frac * sig_idx.size))
        core = sig_idx[n_drop: sig_idx.size - n_drop]
        if core.size == 0:
            core = sig_idx
        for j in range(nV):
            if norm_iso_for[j] not in sig_order:
                continue
            nc = sig_order.index(norm_iso_for[j])
            x = S[core, j] - bmean[j]
            y = S[core, nc] - bmean[nc]
            y_mean = nanmean(y)
            with np.errstate(invalid='ignore'):
                keep = np.isfinite(x) & np.isfinite(y) & (y >= min_den_rel * y_mean)
            n = int(keep.sum())
            if n < 2:
                continue
            r = x[keep] / y[keep]
            mu = nanmean(r)
            sd = nanstd(r)
            with np.errstate(divide='ignore', invalid='ignore'):
                val = 100 * (sd / np.sqrt(n)) / mu if stat == 'sem' else 100 * sd / mu
            if np.isfinite(val):
                sem[ii, j] = val
    return sem


# ----------------------------------------------------------------------------
# Calibration
# ----------------------------------------------------------------------------
def std_name_column(std_vals: pd.DataFrame) -> str:
    cols = list(std_vals.columns)
    for c in cols:
        if str(c).lower() == 'var1':
            return c
    return cols[0]


def std_row_index(std_vals: pd.DataFrame, key: str, set_name: str) -> Optional[int]:
    """Row of the standards table for a standard key (BHV/BCR/BIR need the set
    tag, e.g. 'BHV-GeoRem'; the others match their row key, e.g. GSD1G)."""
    names = std_vals[std_name_column(std_vals)].astype(str).str.upper()
    if key in PRIMARY_STDS:
        tag = ('-' + set_name).upper()
        hit = np.flatnonzero(names.str.contains(key, regex=False) & names.str.contains(tag, regex=False))
    else:
        rk = _norm_str(STANDARD_ROW_KEY.get(key, key))
        nn = [_norm_str(x) for x in names]
        hit = np.flatnonzero([n == rk or n.startswith(rk) for n in nn])
    return int(hit[0]) if hit.size else None


def make_targets_by_set(std_vals: pd.DataFrame, sig_list: Sequence[str], norm_col: str,
                        set_name: str, keys: Sequence[str] = PRIMARY_STDS) -> np.ndarray:
    """makeTargetsBySet generalised -> len(keys) x nSig array of normalised targets
    (NaN where the standard or the element has no value)."""
    if norm_col not in std_vals.columns:
        raise ValueError(f'Normalization column "{norm_col}" not in standards sheet.')
    rows = []
    for key in keys:
        r = std_row_index(std_vals, key, set_name)
        if r is None and key in PRIMARY_STDS:
            raise ValueError(f'Could not find BHV/BCR/BIR rows for set "{set_name}".')
        rows.append(r)
    out = np.full((len(keys), len(sig_list)), np.nan)
    for i, r in enumerate(rows):
        if r is None:
            continue
        den = float(pd.to_numeric(pd.Series([std_vals[norm_col].iloc[r]]), errors='coerce').iloc[0])
        for j, col in enumerate(sig_list):
            if col in std_vals.columns:
                num = float(pd.to_numeric(pd.Series([std_vals[col].iloc[r]]), errors='coerce').iloc[0])
                out[i, j] = num / den
    return out


def linfit_with_intercept(X, Y):
    X = np.asarray(X, float).ravel()
    Y = np.asarray(Y, float).ravel()
    m = np.isfinite(X) & np.isfinite(Y)
    X, Y = X[m], Y[m]
    if X.size < 2:
        return np.nan, np.nan, np.nan
    A = np.column_stack([X, np.ones_like(X)])
    coeff, *_ = np.linalg.lstsq(A, Y, rcond=None)
    yhat = A @ coeff
    sse = np.sum((Y - yhat) ** 2)
    sst = np.sum((Y - Y.mean()) ** 2)
    r2 = np.nan if sst <= EPS else 1 - sse / sst
    return float(coeff[0]), float(coeff[1]), float(r2)


def linfit_through_origin(X, Y):
    X = np.asarray(X, float).ravel()
    Y = np.asarray(Y, float).ravel()
    m = np.isfinite(X) & np.isfinite(Y)
    X, Y = X[m], Y[m]
    if X.size < 2:
        return np.nan, np.nan
    slope = np.sum(X * Y) / max(np.sum(X ** 2), EPS)
    sse = np.sum((Y - slope * X) ** 2)
    sst0 = np.sum(Y ** 2)
    r2 = np.nan if sst0 <= EPS else 1 - sse / sst0
    return float(slope), float(r2)


def weighted_log_slope(X, Y, err_pct, scheme: str = 'ivw'):
    """Sensitivity slope as a weighted mean of log(target/measured) over the
    calibrant points.  'ivw': weights 1/(err^2 + floor^2); 'logmean': equal
    weights.  Returns (slope, R^2 of the through-origin line, n used)."""
    X = np.asarray(X, float).ravel(); Y = np.asarray(Y, float).ravel(); E = np.asarray(err_pct, float).ravel()
    m = np.isfinite(X) & np.isfinite(Y) & (X > 0) & (Y > 0)
    X, Y, E = X[m], Y[m], E[m]
    if X.size < 1:
        return np.nan, np.nan, 0
    w = _weights(E) if scheme == 'ivw' else np.ones_like(X)
    if w.sum() <= 0:
        w = np.ones_like(X)
    slope = float(np.exp(np.sum(w * np.log(Y / X)) / np.sum(w)))
    sse = np.sum((Y - slope * X) ** 2); sst0 = np.sum(Y ** 2)
    r2 = np.nan if sst0 <= EPS else 1 - sse / sst0
    return slope, float(r2), int(X.size)


def parse_cal_exclusions(value) -> List[str]:
    """getCalExclNames on a saved CalExclusions cell: names split on
    whitespace, ignoring blanks and the 'None' placeholder."""
    if is_blank(value):
        return []
    s = str(value).strip()
    if s.lower() == 'none':
        return []
    return [p for p in s.split() if p]


def parse_calibrants(value) -> List[str]:
    """CalibrationSelections 'Calibrants': space/comma list of standard keys."""
    if is_blank(value):
        return list(PRIMARY_STDS)
    toks = [t.strip().upper() for t in re.split(r'[\s,;+]+', str(value)) if t.strip()]
    keys = [t for t in toks if t in STANDARD_TAGS]
    return keys or list(PRIMARY_STDS)


def cal_weighting(cal_prefs: pd.DataFrame, analyte: str, default: str = 'ivw') -> str:
    if analyte in cal_prefs.index and 'Weighting' in cal_prefs.columns and not is_blank(cal_prefs.loc[analyte, 'Weighting']):
        w = str(cal_prefs.loc[analyte, 'Weighting']).lower()
        if w in VALID_WEIGHTINGS:
            return w
    return default


@dataclass
class CalElement:
    analyte: str
    norm_key: str
    set_name: str
    force_zero: bool
    x_plot: np.ndarray
    y_plot: np.ndarray
    names_plot: List[str]
    excluded: List[str]
    slope: float
    intercept: float
    r2: float
    weighting: str = 'ols0'
    calibrants: List[str] = field(default_factory=lambda: list(PRIMARY_STDS))
    err_plot: np.ndarray = field(default_factory=lambda: np.array([]))   # counting error % per point
    auto_excluded: List[str] = field(default_factory=list)              # dropped by the counting-error cut


def cal_norm_key(cal_prefs: pd.DataFrame, analyte: str, norm_map, ui_norm_key: str) -> str:
    key = ''
    if analyte in cal_prefs.index and not is_blank(cal_prefs.loc[analyte, 'NormElement']):
        key = str(cal_prefs.loc[analyte, 'NormElement'])
    if key not in norm_map:
        key = ui_norm_key
    return key


def cal_set_name(cal_prefs: pd.DataFrame, analyte: str, ui_harvard: bool) -> str:
    if analyte in cal_prefs.index:
        v = cal_prefs.loc[analyte, 'StandardSet']
        if isinstance(v, (bool, np.bool_)):
            return 'Harvard' if v else 'GeoRem'
        if is_blank(v):
            return 'GeoRem'
        return 'Harvard' if str(v).lower() == 'harvard' else 'GeoRem'
    return 'Harvard' if ui_harvard else 'GeoRem'


def cal_force_zero(cal_prefs: pd.DataFrame, analyte: str, ui_force_zero: bool) -> bool:
    if analyte in cal_prefs.index:
        v = cal_prefs.loc[analyte, 'ForceInterceptZero']
        if isinstance(v, (bool, np.bool_)):
            return bool(v)
        if isinstance(v, (int, float, np.integer, np.floating)) and not is_blank(v):
            return bool(v)
        if isinstance(v, str) and v.strip().lower() in ('true', 'false', '1', '0'):
            return v.strip().lower() in ('true', '1')
    return ui_force_zero


def norm_index_for(sig_vars: Sequence[str], cal_prefs: pd.DataFrame, norm_map, ui_norm_key: str) -> List[int]:
    """Column index of each analyte's normaliser isotope (the analyte itself when absent)."""
    sig_vars = list(sig_vars)
    out = []
    for j, analyte in enumerate(sig_vars):
        key = cal_norm_key(cal_prefs, analyte, norm_map, ui_norm_key)
        iso = norm_map[key]['iso']
        out.append(sig_vars.index(iso) if iso in sig_vars else j)
    return out


def calibrate_all(lt_corr: np.ndarray, sig_vars: Sequence[str], display_names: Sequence[str],
                  sets: Dict[str, np.ndarray], std_vals: pd.DataFrame, cal_prefs: pd.DataFrame,
                  norm_map, intervals: pd.DataFrame,
                  ui_norm_key: str = 'Ca', ui_force_zero: bool = True, ui_harvard: bool = False,
                  in_memory_excl: Optional[Dict[str, List[str]]] = None,
                  av_total: Optional[np.ndarray] = None, drift_grid: Optional[np.ndarray] = None,
                  cnt_err: Optional[np.ndarray] = None, ratio_drift: Optional[Sequence[bool]] = None,
                  ui_weighting: str = 'ivw') -> Tuple[List[CalElement], np.ndarray]:
    """calibrateAll -> (per-analyte fit info, FinalResult wt%).

    Measured ratio for analyte j: when ``ratio_drift[j]`` (pooled drift) the
    drift factor applies to the ratio, ``(av_total[:, j] / av_total[:, ni]) /
    drift_grid[:, j]``; otherwise the legacy ``lt_corr[:, j] / lt_corr[:, ni]``.
    Weighting per analyte from CalibrationSelections 'Weighting' (ivw, logmean,
    ols0, ols); 'Calibrants' lists the standards used (default BHV BCR BIR).
    Calibrant points with a counting error above MAX_CAL_ERR_PCT are dropped.
    """
    sig_vars = list(sig_vars)
    display_names = list(display_names)
    for k in PRIMARY_STDS:
        if sets[k].size == 0:
            raise ValueError('Need BHV, BCR, and BIR present to calibrate.')
    in_memory_excl = in_memory_excl or {}
    n_samples = lt_corr.shape[0]
    if ratio_drift is None:
        ratio_drift = [False] * len(sig_vars)
    elements: List[CalElement] = []
    meas_all_cols = []

    for j, analyte in enumerate(sig_vars):
        key = cal_norm_key(cal_prefs, analyte, norm_map, ui_norm_key)
        iso = norm_map[key]['iso']
        if iso not in sig_vars:
            raise ValueError(f'Normalization isotope "{iso}" not found in data.')
        ni = sig_vars.index(iso)
        with np.errstate(divide='ignore', invalid='ignore'):
            if ratio_drift[j] and av_total is not None and drift_grid is not None:
                meas_norm = (av_total[:, j] / np.fmax(av_total[:, ni], EPS)) / drift_grid[:, j]
            else:
                meas_norm = lt_corr[:, j] / np.fmax(lt_corr[:, ni], EPS)
        meas_all_cols.append(meas_norm)
        err = ratio_error_pct(cnt_err, j, ni) if cnt_err is not None else np.zeros(n_samples)
        set_name = cal_set_name(cal_prefs, analyte, ui_harvard)
        calibrants = parse_calibrants(cal_prefs.loc[analyte, 'Calibrants']) if (analyte in cal_prefs.index and 'Calibrants' in cal_prefs.columns) else list(PRIMARY_STDS)
        calibrants = [k for k in calibrants if sets.get(k, np.array([], int)).size > 0] or list(PRIMARY_STDS)
        targets = make_targets_by_set(std_vals, sig_vars, iso, set_name, calibrants)
        idx_all = np.concatenate([sets[k] for k in calibrants])
        X = np.concatenate([meas_norm[sets[k]] for k in calibrants])
        Y = np.concatenate([np.repeat(targets[i, j], sets[k].size) for i, k in enumerate(calibrants)])
        E = np.concatenate([err[sets[k]] for k in calibrants])
        names_all = [display_names[i] for i in idx_all]
        fin = np.isfinite(X) & np.isfinite(Y)
        x_plot, y_plot, e_plot = X[fin], Y[fin], E[fin]
        names_plot = [n for n, f in zip(names_all, fin) if f]
        if analyte in in_memory_excl:
            ex = [e for e in in_memory_excl[analyte] if e and e != 'None']
        else:
            ex = parse_cal_exclusions(cal_prefs.loc[analyte, 'CalExclusions']) if analyte in cal_prefs.index else []
        weighting = cal_weighting(cal_prefs, analyte, ui_weighting)
        auto_ex = [n for n, e_ in zip(names_plot, e_plot) if not (e_ <= MAX_CAL_ERR_PCT)] if weighting in ('ivw', 'logmean') else []
        inc = np.array([(n not in ex) and (n not in auto_ex) for n in names_plot], bool)
        if inc.sum() == 0:                     # never drop everything
            inc = np.array([n not in ex for n in names_plot], bool); auto_ex = []
        fz = cal_force_zero(cal_prefs, analyte, ui_force_zero)
        if weighting in ('ivw', 'logmean'):
            slope, r2, _ = weighted_log_slope(x_plot[inc], y_plot[inc], e_plot[inc], weighting)
            intercept = 0.0
            fz = True
        elif weighting == 'ols' or not fz:
            slope, intercept, r2 = linfit_with_intercept(x_plot[inc], y_plot[inc])
            fz = False
        else:
            slope, r2 = linfit_through_origin(x_plot[inc], y_plot[inc])
            intercept = 0.0
        elements.append(CalElement(analyte, key, set_name, fz, x_plot, y_plot, names_plot, ex,
                                   slope, intercept, r2, weighting, calibrants, e_plot, auto_ex))

    final = np.zeros_like(lt_corr, dtype=float)
    for j, el in enumerate(elements):
        oxide_col = norm_map[el.norm_key]['oxideCol']
        if oxide_col in intervals.columns:
            ox = pd.to_numeric(intervals[oxide_col], errors='coerce').to_numpy(float)
        else:
            ox = np.zeros(n_samples)
        ox = np.where(np.isnan(ox), 0.0, ox)
        final[:, j] = (meas_all_cols[j] * el.slope + el.intercept) * ox
    return elements, final


# ----------------------------------------------------------------------------
# Consistency of the standards and down-hole fractionation
# ----------------------------------------------------------------------------
def standard_consistency(elements: List[CalElement], sig_vars: Sequence[str], display_names: Sequence[str],
                         sets: Dict[str, np.ndarray], std_vals: pd.DataFrame, meas_ratio: np.ndarray,
                         cnt_err: Optional[np.ndarray], norm_idx: Sequence[int]) -> pd.DataFrame:
    """Apparent sensitivity of every standard with reference values, relative
    to the calibration slope (1.00 = agrees with the calibrants). Columns per
    standard: ratio and its counting error %; rows = analytes."""
    rows = {}
    keys = [k for k in STANDARD_TAGS if sets.get(k, np.array([], int)).size > 0]
    for j, el in enumerate(elements):
        r = {}
        for k in keys:
            ri = std_row_index(std_vals, k, el.set_name)
            if ri is None or el.analyte not in std_vals.columns:
                continue
            norm_iso = sig_vars[norm_idx[j]]
            tgt = pd.to_numeric(pd.Series([std_vals[el.analyte].iloc[ri], std_vals[norm_iso].iloc[ri]]), errors='coerce').to_numpy(float)
            if not np.all(np.isfinite(tgt)) or tgt[1] == 0 or tgt[0] <= 0:
                continue
            idx = sets[k]
            x = meas_ratio[idx, j]
            e = ratio_error_pct(cnt_err, j, norm_idx[j])[idx] if cnt_err is not None else np.zeros(idx.size)
            w = _weights(e); good = np.isfinite(x) & (x > 0) & (w > 0)
            if not good.any():
                continue
            mu = np.sum(w[good] * x[good]) / np.sum(w[good])
            r[(k, 'ratio')] = mu * el.slope / (tgt[0] / tgt[1])
            r[(k, 'err%')] = 1.0 / np.sqrt(np.sum(w[good])) if np.isfinite(mu) else np.nan
        rows[el.analyte] = r
    return pd.DataFrame(rows).T


def downhole_slopes(windows: List[pd.DataFrame], intervals: pd.DataFrame, sig_order: Sequence[str],
                    norm_idx: Sequence[int], idx: np.ndarray) -> pd.DataFrame:
    """Slope of ln(analyte/normaliser) vs time inside the signal window, % per
    10 s, per analysis in ``idx``; returns the median and the spread over them."""
    sig_order = list(sig_order)
    out = {}
    per = []
    for ii in idx:
        T = windows[ii]
        t = T['Time'].to_numpy(float)
        S = T[sig_order].to_numpy(float)
        s0, s1 = _iv(intervals, ii, 'signal_start'), _iv(intervals, ii, 'signal_stop')
        b1, e1 = _iv(intervals, ii, 'back_start1'), _iv(intervals, ii, 'back_stop1')
        with np.errstate(invalid='ignore'):
            mb = (t >= b1) & (t < e1) if np.isfinite(b1) and np.isfinite(e1) else np.zeros_like(t, bool)
            ms = (t >= s0) & (t < s1)
        if ms.sum() < 6:
            continue
        bk = nanmean_rows(S[mb]) if mb.any() else np.zeros(S.shape[1])
        net = S[ms] - bk
        tt = t[ms] - t[ms].mean()
        row = np.full(len(sig_order), np.nan)
        for j in range(len(sig_order)):
            ni = norm_idx[j]
            if ni == j:
                row[j] = 0.0
                continue
            with np.errstate(divide='ignore', invalid='ignore'):
                y = np.log(net[:, j] / net[:, ni])
            ok = np.isfinite(y)
            if ok.sum() >= 6:
                row[j] = 1000.0 * np.polyfit(tt[ok], y[ok], 1)[0]   # %/10 s
        per.append(row)
    if not per:
        return pd.DataFrame(index=sig_order, columns=['slope_pct_per_10s', 'spread_pct_per_10s', 'n'])
    P = np.array(per)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', category=RuntimeWarning)
        med = np.nanmedian(P, axis=0)
        mad = 1.4826 * np.nanmedian(np.abs(P - med), axis=0)
    return pd.DataFrame({'slope_pct_per_10s': med, 'spread_pct_per_10s': mad, 'n': np.isfinite(P).sum(axis=0)}, index=sig_order)


# ----------------------------------------------------------------------------
# Secondary standards summary (right-hand table of the Calibration tab)
# ----------------------------------------------------------------------------
def _norm_str(s) -> str:
    return re.sub(r'[^A-Z0-9]', '', str(s).upper())


def secondary_summary(display_names: Sequence[str], analyte: str, sig_vars: Sequence[str],
                      final: np.ndarray, std_vals: Optional[pd.DataFrame],
                      set_name: str) -> List[Tuple]:
    """buildSecondarySummaryForAnalyte -> rows (Std, N, Mean, %RSD, Target, Mean/Target, %Bias)."""
    if std_vals is None or len(std_vals) == 0:
        return []
    name_col = std_name_column(std_vals)
    names_norm = [_norm_str(s) for s in std_vals[name_col].astype(str)]
    set_norm = _norm_str(set_name)
    primaries = ('BHV', 'BCR', 'BIR')

    groups: Dict[str, List[int]] = {}
    for i, nm in enumerate(display_names):
        low = str(nm).lower()
        if 'bhv' in low:
            key = 'BHV'
        elif 'bcr' in low:
            key = 'BCR'
        elif 'bir' in low:
            key = 'BIR'
        elif 've32' in low:
            key = 'VE32'
        elif 'otherst' in low:
            key = 'OtherST'
        elif 'sth' in low:
            key = 'STHS'
        elif 'gsd1g' in low or 'gsd' in low:
            key = 'GSD1G'
        elif 'gse1g' in low or 'gse' in low:
            key = 'GSE1G'
        elif 'gor_128' in low or 'gor-128' in low or 'gor 128' in low:
            key = 'GOR-128'
        else:
            key = str(nm)
        kn = _norm_str(key)
        if key in primaries:
            mask = [kn in n and set_norm in n for n in names_norm]
        else:
            mask = [n == kn or kn in n for n in names_norm]
        if not any(mask):
            continue
        groups.setdefault(key, []).append(i)

    sig_vars = list(sig_vars)
    if analyte not in sig_vars:
        return []
    j = sig_vars.index(analyte)
    rows = []
    for key in sorted(groups):            # containers.Map keys come back sorted
        idxs = groups[key]
        vec = final[idxs, j]
        mu, sd = nanmean(vec), nanstd(vec)
        n = int(np.isfinite(vec).sum())
        rsd = 100 * sd / mu if (np.isfinite(mu) and mu != 0 and n > 1) else np.nan
        tgt = ratio = bias = np.nan
        if analyte in std_vals.columns:
            sn = _norm_str(key)
            if key in primaries:
                hits = [k for k, n_ in enumerate(names_norm) if sn in n_ and set_norm in n_]
            else:
                hits = [k for k, n_ in enumerate(names_norm) if n_ == sn or sn in n_]
            if hits:
                v = pd.to_numeric(pd.Series([std_vals[analyte].iloc[hits[0]]]), errors='coerce').iloc[0]
                tgt = float(v) if pd.notna(v) else np.nan
        if np.isfinite(mu) and np.isfinite(tgt) and tgt != 0:
            ratio = mu / tgt
            bias = 100 * (mu - tgt) / tgt
        rows.append((key, n, mu, rsd, tgt, ratio, bias))
    return rows
