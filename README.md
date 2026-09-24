# LA_Quickcheck

Portable first-pass data checker for LA-ICP-MS runs (WHOI/UMass reduction scheme).
Point it at a run, and it picks reduction parameters from the data, reduces it,
and produces calibration-curve and standard-recovery plots with an HTML report.
Runs from a thumb drive on a Windows PC with no Python installed.

This branch is the **precision fork**: the same file layouts and the same
MATLAB-compatible numerics as before (selectable per element), plus a pooled,
count-weighted drift model, an inverse-variance calibration, longer signal
windows, counting-statistics floors in the report, and standard-consistency /
down-hole-fractionation tables. See "What changed in the precision fork" below.

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
   DwellTimes.xlsx              dwell per isotope read from the Qtegra files (counting statistics)
   DriftSelections.xlsx         drift standard (ALL = pooled), model, exclusions per analyte
   CalibrationSelections.xlsx   normaliser, intercept, standard set, exclusions, weighting, calibrants per analyte
   DriftCorrected_All.xlsx      drift-corrected background-subtracted signals
   ReducedDataExport.xlsx       FinalResult, Signal_SEM_Norm, Counting_Err_Pct, snapshots of all settings
   StandardConsistency.xlsx     apparent sensitivity of every standard relative to the calibration
   DownholeSlopes.xlsx          analyte/normaliser slope inside the signal window, % per 10 s
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
| Background / signal | From the 25Mg trace of each analysis: background 1 s → laser-on − 2 s; signal laser-on + 2 s for (standards' median ablation − 4 s), the same length for every sample. The first second carries a surface Pb/Cu spike and a sweep-timing artefact; the primaries are counting-limited, so the window is the whole ablation. An unknown that ablated for less gets its own (ablation − 4 s) window and is flagged in the report. Second background only when ≥ 12 s of clean blank follow the ablation. |
| Drift | `ALL` / `auto`: every recognised standard in every bracket is pooled. For each element the ratio to its normaliser is divided by that standard's own count-weighted run mean, then averaged per bracket with weights 1/(counting error² + 0.5 %²). pchip through the bracket means when they are good to 1.5 %, quadratic otherwise. The factor applies to the ratio, so the normaliser isotope itself carries no correction. Set a single key (BHV …) and a MATLAB method per element in `DriftSelections.xlsx` to get the legacy single-standard model on the raw signal. |
| Normaliser | Chosen per analyte from the run's own standards: every available normaliser (Al, Ca, Si) is tried on the whole run and the one with the lowest score, replicate RSD and cross-standard spread of apparent sensitivity in quadrature, is kept (`data/NormaliserChoice.xlsx`); the Al default survives unless another beats it by 0.2 %. On `09_21_26_50um`: 22 analytes on Al, 7 on Ca, 5 on Si (Na, Ni, Cu, Ce, Pb). |
| Calibration | Slope = inverse-variance weighted mean of log(target / measured) over the BHV/BCR/BIR replicates (`Weighting` = `ivw`; `logmean`, `ols0` (legacy through-zero OLS), `ols` also available per element). Replicates with > 10 % counting error are dropped from the fit. `Calibrants` lists the standards used (any of BHV BCR BIR GSD GSE STH); default BHV BCR BIR except Cu, which uses BHV BIR because the BCR-2G Cu value is inconsistent with every other standard (`core.DEFAULT_CALIBRANTS`). GeoRem values. |
| Exclusions | A calibrant replicate more than 4 robust σ and more than 10 % from its own standard's other replicates is excluded (never more than a third of a standard's points). Standard-to-standard offsets are reported in the consistency table, not "fixed". |
| Oxide wt% | Autofilled for BHV, BCR, BIR, GSD, GSE, StHS, VE32, GOR-128. Unknowns read 0 until you enter their Al2O3 / CaO in `data/Intervals.xlsx` and re-run. |

Standards are recognised by name: `bhv`, `bcr`, `bir`, `gsd`, `gse`, `sth`, `ve32`, `gor_128`
(case-insensitive substring). Reference values come from `LAICPMS_Standard_Values.xlsx`
in this repository, which is bundled into the executable; a copy is placed in every
`data/` folder for the record.

**GSE-1G values are provisional.** The `GSE1G-provisional` row was derived from run
`09_21_26_50um` calibrated against BHV/BCR/BIR with this fork (oxides set to the GSD-1G
matrix). It makes GSE-1G count for drift and lets it appear in the plots; its bias column is
self-referential until the row is replaced by GeoRem preferred values (Jochum et al. 2005).
Do not list GSE as a calibrant before that.

## What changed in the precision fork

Measured on `09_21_26_50um` (six standards × 7 brackets, 34 isotopes at 10 ms dwell).
Between-replicate RSD, median over the trace elements that are well measured in every
standard (the same code path an unknown sees; GSD/GSE are inside the pooled drift model,
so their numbers are slightly in-sample):

| | BHV | BCR | BIR | GSD | GSE | StHs |
|---|---|---|---|---|---|---|
| previous defaults (BHV-only quadratic, 1 + 20 s, OLS through zero) | 1.3 | 2.4 | 3.2 | 1.25 | 1.36 | 3.0 |
| this fork (pooled weighted pchip, 2 + 28 s, inverse-variance) | 1.3 | 1.5 | 3.1 | 0.56 | 0.39 | 2.3 |

Median |bias| of the same elements: BHV 0.75 %, BCR 0.9 %, BIR 1.6 %, GSD 1.1 %, StHs 3.9 %.
The primaries sit on their counting-statistics floor (10 ms × 34 masses = 2.6 % duty cycle
per isotope; 175Lu in BHVO-2G is ~5 % per analysis by counting alone), so further gains
there need longer dwells on the low-abundance isotopes, not reduction changes.

What the report adds: the counting floor next to every RSD (bias hidden where the floor
is > 5 %), the consistency of every standard against the calibration (a standard off for
one element = that reference value; off for every element = its internal-standard oxide or
normaliser), the standards' down-hole fractionation slopes, and the list of unknowns that
ablated for less than the standards.

Two reference-value checks done on this run: the GeoRem rows beat the Harvard rows (median
spread among the primaries 3.7 % vs 6.4 %, driven by the Harvard BIR-1G row; Harvard is
better only for Sr and Pb, selectable per element via `StandardSet`). For Pb, clipping more
of the early signal does not help: the surface spike is confined to the first second, and the
±9 % BHV/BCR disagreement is the BCR-2G reference value (calibrating Pb on `BHV BIR` brings
GSD to −1 % and StHs to +8 %).

Down-hole fractionation was assessed on the standards and nothing is corrected. For
full-length analyses no scheme helps: the window mean beats a linear-fit intercept two-fold in
precision, and the replicate scatter is uncorrelated with each spot's decay, slope or yield
(all 42 standard ablations decay to 0.62 ± 0.02 of their initial Al). Referring short-window
analyses to the standards' window with the standards' slopes did halve the bias in a
leave-one-standard-out test (10 s windows: 1.1 → 0.6 %), but it assumes an unknown fractionates
like the glass standards, which is not safe for inclusions, so it was not adopted. Analyses whose
window is shorter than the standards' are flagged in the report with the slopes in
`data/DownholeSlopes.xlsx` for the user to judge.

Legacy behaviour is one setting away: `Standard` = a single key and a MATLAB method in
`DriftSelections.xlsx`, `Weighting` = `ols0` in `CalibrationSelections.xlsx`, and the
window of your choice in `Intervals.xlsx` reproduce the previous outputs exactly
(checked to 0 relative difference on this run).

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
   core.py                  numerical steps (windows, background, counting errors, pooled + legacy drift,
                            weighted + legacy calibration, consistency, down-hole slopes, SEM)
   session.py               orchestration (load -> drift -> calibrate -> export)
   io_xlsx.py               all XLSX reading / writing, MATLAB-compatible layouts
   interp.py, names.py      MATLAB-compatible interpolation and naming rules
plot_calibrations.py        calibration-curve plots
plot_standard_recovery.py   reference vs reduced value plots with ±2 SEM error bars
build_datachecker.py        PyInstaller build
LAICPMS_Standard_Values.xlsx  reference values (GeoRem / Harvard rows for BHV, BCR, BIR; secondaries)
```

The numerical core is a transcription of the MATLAB `LaserCalTool_Beta_17.m` and was
verified to reproduce that tool's exports to ~1e-13 relative on three reference runs; the
legacy code paths are unchanged and selectable per element (see above).

Tests: `python -m pytest tests` (pooled drift, counting errors, weighted calibration).
