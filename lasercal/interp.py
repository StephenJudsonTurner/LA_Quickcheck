"""MATLAB-compatible 1-D interpolation and polynomial fitting.

The drift tab of the MATLAB tool calls ``interp1(xg, yg, tAll, method, 'extrap')``
with method in {linear, pchip, makima, spline} and ``polyfit``/``polyval`` for
poly2 / poly3.  SciPy implements the same algorithms (PCHIP uses the same
end-point and Fritsch-Butland slope rules as MATLAB, CubicSpline defaults to
not-a-knot, Akima1DInterpolator(method='makima') is MATLAB's modified Akima),
so results agree to round-off.  This module wraps them so the edge cases
(two points -> straight line, extrapolation) behave like MATLAB.
"""
from __future__ import annotations

import numpy as np
from scipy.interpolate import PchipInterpolator, CubicSpline, Akima1DInterpolator

METHODS = ("none", "linear", "pchip", "makima", "spline", "poly2", "poly3")


def _linear_extrap(xg, yg, xq):
    """interp1(...,'linear','extrap'): straight-line extrapolation from end segments."""
    xg = np.asarray(xg, float)
    yg = np.asarray(yg, float)
    xq = np.asarray(xq, float)
    if xg.size < 2:
        raise ValueError("interp1 requires at least two sample points.")
    out = np.interp(xq, xg, yg)
    lo = xq < xg[0]
    hi = xq > xg[-1]
    if lo.any():
        s = (yg[1] - yg[0]) / (xg[1] - xg[0])
        out[lo] = yg[0] + s * (xq[lo] - xg[0])
    if hi.any():
        s = (yg[-1] - yg[-2]) / (xg[-1] - xg[-2])
        out[hi] = yg[-1] + s * (xq[hi] - xg[-1])
    return out


def interp1_extrap(xg, yg, xq, method: str):
    """Equivalent of MATLAB ``interp1(xg, yg, xq, method, 'extrap')``.

    ``xg`` must be strictly increasing (the MATLAB tool always passes
    group-mean start times, which are).  Raises ValueError for < 2 points,
    which is what MATLAB does (the caller catches and falls back).
    """
    xg = np.asarray(xg, float).ravel()
    yg = np.asarray(yg, float).ravel()
    xq = np.asarray(xq, float).ravel()
    method = method.lower()
    if xg.size < 2:
        raise ValueError("interp1 requires at least two sample points.")
    if np.any(np.diff(xg) <= 0):
        order = np.argsort(xg)
        xg, yg = xg[order], yg[order]
        if np.any(np.diff(xg) <= 0):
            raise ValueError("interp1 sample points must be unique.")
    if method == "linear" or xg.size == 2:
        # MATLAB: pchip/spline/makima with exactly two points reduce to a line.
        return _linear_extrap(xg, yg, xq)
    if method == "pchip":
        f = PchipInterpolator(xg, yg, extrapolate=True)
    elif method == "spline":
        f = CubicSpline(xg, yg, bc_type="not-a-knot", extrapolate=True)
    elif method == "makima":
        f = Akima1DInterpolator(xg, yg, method="makima", extrapolate=True)
    else:
        raise ValueError(f"Unknown interpolation method '{method}'")
    return f(xq)


def polyfit_polyval(xg, yg, xq, deg: int):
    """MATLAB ``polyval(polyfit(xg, yg, deg), xq)`` (least squares on the
    unscaled Vandermonde matrix, as MATLAB does)."""
    xg = np.asarray(xg, float).ravel()
    yg = np.asarray(yg, float).ravel()
    xq = np.asarray(xq, float).ravel()
    V = np.vander(xg, deg + 1)
    p, *_ = np.linalg.lstsq(V, yg, rcond=None)
    return np.polyval(p, xq)
