# LA_Quickcheck

Portable first-pass data checker for LA-ICP-MS runs (WHOI/UMass reduction scheme).
Point it at a run, and it picks reduction parameters from the data, reduces it,
and produces calibration-curve and standard-recovery plots with an HTML report.
Runs from a thumb drive on a Windows PC with no Python installed.

## Download and run (no Python needed)

1. Get `DataChecker_win64.zip` from the **Releases** page and unzip it anywhere
   (thumb drive is fine). Keep `DataChecker.exe` and `_internal/` together.
2. Double-click `DataChecker.exe`.
3. Browse to a folder of Thermo Qtegra per-analysis CSV files (one per spot),
   or to a data folder holding one data workbook (`Time` + isotope columns) and
   `RunOrder.xlsx`.
4. Press **Run**, then **Open report**.

Everything is written to `<input>_checked/` next to the input:

```
REPORT.html                     what was chosen, per-standard bias / RSD table, contact sheets
data/                           all settings and results as XLSX (editable)
   <run>_data.xlsx              concatenated traces (CSV input) or a copy of the workbook
   RunOrder.xlsx, SampleWindows.xlsx, SourceFiles.xlsx
   Intervals.xlsx               background / signal windows + internal-standard oxide wt%
   DriftSelections.xlsx         drift standard, model, exclusions per analyte
   CalibrationSelections.xlsx   normaliser, intercept, standard set, exclusions per analyte
   DriftCorrected_All.xlsx      drift-corrected background-subtracted signals
   ReducedDataExport.xlsx       FinalResult, Signal_SEM_Norm, snapshots of all settings
plots/calibration_curves/       one PNG per analyte + index.png + calibration_stats.xlsx
plots/standard_recovery/        one PNG per analyte + index.png + recovery_stats.xlsx
```

To refine a run, edit the XLSX files in `data/` (same layouts as the MATLAB
LaserCalTool), choose **Re-run a *_checked folder** in the app and press Run.

Command line: `DataChecker.exe <input folder>` or `DataChecker.exe --rerun <..._checked>`.

## What it decides automatically

| Choice | Rule |
|---|---|
| Sample windows | CSV input: each file is one analysis. Workbook input: the MATLAB threshold search (25Mg ≥ 1000, 80-row gap). |
| Background / signal | From the 25Mg trace of each analysis: background 1 s → laser-on − 2 s; signal laser-on + 1 s for min(20 s, ablation − 5 s), same length for every sample; second background only when ≥ 12 s of clean blank follow the ablation. (A window sweep on a 50 µm run showed precision flattening beyond ~20 s and no sensitivity to the delay.) |
| Drift | Standard: BHV if present, else the most frequent of BCR/GSD/BIR. Model: pchip through bracket means when the standard was run in consecutive pairs; quadratic when run singly in ≥ 4 brackets (an interpolant would pin every point to its mean); linear for 2–3 brackets; none for 1. |
| Calibration | BHV, BCR, BIR against GeoRem values, normalised to Al (Ca for the Al isotope), intercept forced through zero. |
| Exclusions | A calibrant replicate more than 4 robust σ and more than 10 % from its own standard's other replicates is excluded (never more than a third of a standard's points). Standard-to-standard offsets are reported, not "fixed". |
| Oxide wt% | Autofilled for BHV, BCR, BIR, GSD, StHS, VE32, GOR-128. Unknowns read 0 until you enter their Al2O3 / CaO in `data/Intervals.xlsx` and re-run. |

Standards are recognised by name: `bhv`, `bcr`, `bir`, `gsd`, `sth`, `ve32`, `gor_128`
(case-insensitive substring). Reference values come from `LAICPMS_Standard_Values.xlsx`
in this repository, which is bundled into the executable; a copy is placed in every
`data/` folder for the record.

## Building the executable yourself

Python 3.10 with `numpy`, `pandas`, `scipy ≥ 1.13`, `openpyxl`, `matplotlib`, `pyinstaller`
(see `requirements.txt`):

```bash
python build_datachecker.py
```

produces `dist/DataChecker/` (one folder, about 290 MB). Zip that folder for distribution.

Running from source without building: `python datachecker_app.py [input]`.

## Layout

```
datachecker_app.py          window / command-line entry point
lasercal/                   the reduction library
   checker.py               automatic choices, run, HTML report
   qtegra.py                Qtegra per-analysis CSV -> tool inputs
   core.py                  numerical steps (windows, background, drift, calibration, SEM)
   session.py               orchestration (load -> drift -> calibrate -> export)
   io_xlsx.py               all XLSX reading / writing, MATLAB-compatible layouts
   interp.py, names.py      MATLAB-compatible interpolation and naming rules
plot_calibrations.py        calibration-curve plots
plot_standard_recovery.py   reference vs reduced value plots with ±2 SEM error bars
build_datachecker.py        PyInstaller build
LAICPMS_Standard_Values.xlsx  reference values (GeoRem / Harvard rows for BHV, BCR, BIR; secondaries)
```

The numerical core is a transcription of the MATLAB `LaserCalTool_Beta_17.m` and was
verified to reproduce that tool's exports to ~1e-13 relative on three reference runs.
