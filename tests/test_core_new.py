"""Unit tests for the pooled drift model, counting errors and weighted calibration."""
import numpy as np
import pandas as pd
import pytest

from lasercal import core


def _sets(n_brackets=4, stds=('BHV', 'BCR', 'BIR', 'GSD'), n_unknown=5):
    """Run order: each bracket = one analysis of every standard, then unknowns."""
    names, sets = [], {k: [] for k in stds}
    for b in range(n_brackets):
        for k in stds:
            sets[k].append(len(names)); names.append(f'{k}_{b + 1}')
        for u in range(n_unknown):
            names.append(f'UNK_{b}_{u}')
    return names, {k: np.array(v, int) for k, v in sets.items()}


def test_standard_brackets_split_back_to_back():
    names, sets = _sets(n_brackets=2, n_unknown=0)     # brackets adjacent
    groups = core.standard_brackets(sets)
    assert len(groups) == 2 and groups[0].tolist() == [0, 1, 2, 3]


def test_counting_error_matches_poisson():
    n_sig, n_back, dwell = 50, 60, 0.01
    gross, back = 1000.0, 10.0
    st = core.SignalStats(np.array([[gross - back]]), np.array([[back]]), np.array([[gross]]),
                          np.array([n_sig]), np.array([n_back]))
    e = core.counting_error_pct(st, np.array([dwell]))
    expect = 100 * np.sqrt(gross / (n_sig * dwell) + back / (n_back * dwell)) / (gross - back)
    assert e[0, 0] == pytest.approx(expect)
    st.av_total[0, 0] = -1.0
    assert not np.isfinite(core.counting_error_pct(st, np.array([dwell]))[0, 0])


def test_pooled_drift_recovers_common_drift():
    rng = np.random.default_rng(0)
    names, sets = _sets(n_brackets=5)
    n = len(names)
    t = np.arange(n, dtype=float) * 60.0
    true = 1 + 0.05 * np.sin(t / t.max() * 3)             # common multiplicative drift
    conc = {'BHV': 1.0, 'BCR': 3.0, 'BIR': 0.1, 'GSD': 10.0}
    ratio = np.ones(n); err = np.full(n, 0.3)
    for k, idx in sets.items():
        ratio[idx] = conc[k] * true[idx] * (1 + rng.normal(0, 0.003, idx.size))
    err[sets['BIR']] = 30.0                              # BIR is noise: must get ~no weight
    ratio[sets['BIR']] *= (1 + rng.normal(0, 0.3, 5))
    pooled = core.pooled_bracket_means(ratio, err, t, sets, list(sets), np.array([], int))
    assert pooled is not None and pooled.group_times.size == 5
    curve = core.drift_curve('pchip', pooled.group_times, pooled.group_vals, t)
    # compare shape (each normalised to its mean) on the unknowns
    unk = np.array([i for i in range(n) if 'UNK' in names[i]])
    rel_true = true[unk] / true[unk].mean(); rel_fit = curve[unk] / curve[unk].mean()
    assert np.max(np.abs(rel_fit / rel_true - 1)) < 0.01
    assert np.all(pooled.group_se_pct < 1.0)
    assert pooled.group_se_pct.max() < 1.5    # the BIR noise did not leak in


def test_compute_drift_grid_pooled_vs_single():
    names, sets = _sets(n_brackets=3)
    n = len(names); t = np.arange(n, dtype=float)
    av = np.ones((n, 2)); av[:, 0] = 2.0 * (1 + 0.1 * t / t.max())
    err = np.full((n, 2), 0.2)
    prefs = pd.DataFrame({'Standard': ['ALL', 'ALL'], 'Interpolation': ['auto', 'auto'], 'Exclusions': ['', '']},
                         index=['x88Sr', 'x27Al'])
    g = core.compute_drift_grid(av, t, ['x88Sr', 'x27Al'], prefs, sets, cnt_err=err, norm_idx=[1, 1])
    assert np.allclose(g[:, 1], 1.0)                      # normaliser untouched
    corrected = (av[:, 0] / av[:, 1]) / g[:, 0]
    assert np.std(corrected) / np.mean(corrected) < 0.01  # drift removed
    prefs.loc['x88Sr', 'Standard'] = 'BHV'; prefs.loc['x88Sr', 'Interpolation'] = 'linear'
    g2 = core.compute_drift_grid(av, t, ['x88Sr', 'x27Al'], prefs, sets, cnt_err=err, norm_idx=[1, 1])
    assert not np.allclose(g2[:, 0], g[:, 0])              # legacy path is different code


def test_weighted_log_slope_ignores_noisy_points():
    X = np.array([1.0, 2.0, 0.001]); Y = np.array([2.0, 4.0, 0.01])      # third point 5x off
    E = np.array([0.5, 0.5, 50.0])
    s, r2, n = core.weighted_log_slope(X, Y, E, 'ivw')
    assert s == pytest.approx(2.0, rel=1e-3) and n == 3
    s_eq, _, _ = core.weighted_log_slope(X, Y, E, 'logmean')
    assert s_eq > 2.5


def test_parse_calibrants_and_weighting():
    assert core.parse_calibrants('BHV, BCR GSD') == ['BHV', 'BCR', 'GSD']
    assert core.parse_calibrants(None) == list(core.PRIMARY_STDS)
    assert core.parse_calibrants(None, 'x65Cu') == ['BHV', 'BIR']      # BCR-2G Cu dropped by default
    assert core.parse_calibrants('BHV BCR BIR', 'x65Cu') == ['BHV', 'BCR', 'BIR']  # explicit wins
    assert core.default_calibrants('x63Cu') == ['BHV', 'BIR'] and core.default_calibrants('x88Sr') == list(core.PRIMARY_STDS)
    prefs = pd.DataFrame({'Weighting': ['OLS0', 'bogus']}, index=['a', 'b'])
    assert core.cal_weighting(prefs, 'a') == 'ols0' and core.cal_weighting(prefs, 'b') == 'ivw'


def test_auto_method_switch():
    assert core._auto_method(7, 0.5) == 'pchip'
    assert core._auto_method(7, 3.0) == 'poly2'
    assert core._auto_method(2, 0.1) == 'linear'
    assert core._auto_method(1, 0.1) == 'none'


def test_robust_slope_downweights_bad_reference():
    # 14 points agreeing at slope 2, 7 points (one "standard") 15 % off: robust ivw ignores them, ivw0 does not
    X = np.concatenate([np.ones(7), 2 * np.ones(7), np.ones(7)])
    Y = np.concatenate([2 * np.ones(7), 4 * np.ones(7), 2 * 0.85 * np.ones(7)])
    E = np.full(21, 0.5)
    s_rob, _, _ = core.weighted_log_slope(X, Y, E, 'ivw')
    s_plain, _, _ = core.weighted_log_slope(X, Y, E, 'ivw0')
    assert s_rob == pytest.approx(2.0, rel=1e-3)
    assert 1.85 < s_plain < 1.97


def test_common_offsets_centre_to_zero():
    dev = {'A': {f'e{i}': 0.03 for i in range(6)}, 'B': {f'e{i}': -0.01 for i in range(6)}}
    err = {'A': {f'e{i}': 1.0 for i in range(6)}, 'B': {f'e{i}': 1.0 for i in range(6)}}
    c = core.standard_common_offsets(dev, err)
    assert c['A'] == pytest.approx(0.02) and c['B'] == pytest.approx(-0.02)
    err['A'] = {f'e{i}': 9.0 for i in range(6)}       # too noisy to define an offset
    c = core.standard_common_offsets(dev, err)
    assert c['A'] == pytest.approx(0.005) and c['B'] == pytest.approx(-0.005)
