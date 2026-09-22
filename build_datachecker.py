"""Build the portable DataChecker with PyInstaller.

    C:/Users/temp_/miniconda3/envs/thermoengine/python.exe build_datachecker.py

Produces dist/DataChecker/DataChecker.exe (one-folder
build: copy the whole DataChecker folder to a thumb drive; no Python needed
on the target PC).  The standards table LAICPMS_Standard_Values.xlsx in this folder is bundled.
"""
import os
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
STD = os.path.join(HERE, 'LAICPMS_Standard_Values.xlsx')
DIST = os.path.join(HERE, 'dist')
WORK = os.path.join(os.environ.get('TEMP', HERE), 'pyinstaller_datachecker')
TMP_DIST = os.path.join(WORK, 'dist')   # build outside OneDrive (it locks folders during deletion), then swap in

cmd = [sys.executable, '-m', 'PyInstaller', '--noconfirm', '--clean', '--onedir', '--windowed',
       '--name', 'DataChecker',
       '--distpath', TMP_DIST, '--workpath', WORK, '--specpath', WORK,
       '--add-data', f'{STD}{os.pathsep}.',
       '--hidden-import', 'plot_calibrations', '--hidden-import', 'plot_standard_recovery',
       '--hidden-import', 'lasercal.checker', '--hidden-import', 'lasercal.qtegra',
       '--collect-submodules', 'lasercal',
       '--exclude-module', 'PySide6', '--exclude-module', 'PyQt5', '--exclude-module', 'PyQt6',
       '--exclude-module', 'pytest', '--exclude-module', 'IPython', '--exclude-module', 'notebook',
       '--exclude-module', 'thermoengine', '--exclude-module', 'PIL.ImageQt',
       '--paths', HERE,
       os.path.join(HERE, 'datachecker_app.py')]
print(' '.join(cmd))
subprocess.run(cmd, check=True, cwd=HERE)

# swap the fresh build into place: rename the old folder aside if it cannot be deleted
import time
dest = os.path.join(DIST, 'DataChecker')
os.makedirs(DIST, exist_ok=True)
if os.path.isdir(dest):
    try:
        shutil.rmtree(dest)
    except PermissionError:
        aside = dest + '_old_' + time.strftime('%H%M%S')
        os.rename(dest, aside)
        print('old build could not be deleted (file lock); moved to', aside)
shutil.copytree(os.path.join(TMP_DIST, 'DataChecker'), dest)

# a README next to the exe
readme = os.path.join(DIST, 'DataChecker', 'README.txt')
with open(readme, 'w', encoding='utf-8') as f:
    f.write("""LA-ICP-MS Data Checker (portable)
==================================
Copy this whole folder to the lab PC or run it from the thumb drive.
No Python installation is needed.

1. Double-click DataChecker.exe.
2. Browse to a folder of Qtegra CSV files (one per analysis), or to a data
   folder that holds one data workbook plus RunOrder.xlsx.
3. Press Run. Parameters are chosen automatically (windows from the laser-on
   time, drift model from how BHV was bracketed, Al normalisation, GeoRem,
   gross-outlier exclusion) and everything is written to <input>_checked\\
   next to the input:  REPORT.html, data\\ (all XLSX settings + results), plots\\.
4. Press "Open report".

To refine: edit data\\Intervals.xlsx, DriftSelections.xlsx or
CalibrationSelections.xlsx in the _checked folder (same layouts as the MATLAB
tool), choose "Re-run a *_checked folder" and press Run.

Command line (optional):
  DataChecker.exe <input folder>            new check
  DataChecker.exe --rerun <..._checked>     re-run edited settings
The standards table used is LAICPMS_Standard_Values.xlsx inside _internal\\;
a copy is placed in every data\\ folder for the record.
""")
print('built:', os.path.join(DIST, 'DataChecker', 'DataChecker.exe'))
