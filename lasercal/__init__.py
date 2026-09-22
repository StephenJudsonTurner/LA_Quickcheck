"""
lasercal - Python re-implementation of the WHOI LA-ICP-MS reduction tool
(LaserCalTool_Beta_17.m).

Package layout
--------------
names.py    MATLAB-compatible name mangling (readtable / makeValidName rules,
            sample-name cleaning and uniquifying).
interp.py   MATLAB-compatible interp1 / polyfit wrappers.
core.py     Pure numerical steps: sample-start detection, windowing,
            background-corrected averages, drift grid, calibration fits,
            SEM%, secondary-standard summary.  No file or GUI code.
io_xlsx.py  Reading / writing every XLSX file with the exact layouts the
            MATLAB tool uses.
session.py  ReductionSession: the MATLAB startup sequence + the "Calculate
            drifts", "Compute", "Export" actions, with the same fall-backs
            to UI state that the MATLAB code has (no GUI yet; see README).
"""
__version__ = "1.0.0"
