LA-ICP-MS Data Checker (portable)
==================================
Copy this whole folder to the lab PC or run it from the thumb drive.
No Python installation is needed.

1. Double-click DataChecker.exe.
2. Browse to a folder of Qtegra CSV files (one per analysis), or to a data
   folder that holds one data workbook plus RunOrder.xlsx.
3. Press Run. Parameters are chosen automatically (windows from the laser-on
   time, drift model from how BHV was bracketed, Al normalisation, GeoRem,
   gross-outlier exclusion) and everything is written to <input>_checked\
   next to the input:  REPORT.html, data\ (all XLSX settings + results), plots\.
4. Press "Open report".

To refine: edit data\Intervals.xlsx, DriftSelections.xlsx or
CalibrationSelections.xlsx in the _checked folder (same layouts as the MATLAB
tool), choose "Re-run a *_checked folder" and press Run.

Command line (optional):
  DataChecker.exe <input folder>            new check
  DataChecker.exe --rerun <..._checked>     re-run edited settings
The standards table used is LAICPMS_Standard_Values.xlsx inside _internal\;
a copy is placed in every data\ folder for the record.
