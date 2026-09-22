"""LA-ICP-MS Data Checker - stand-alone entry point.

Double-click (or run with no arguments) for a small window: pick a folder of
Qtegra CSV files or a data folder (workbook + RunOrder.xlsx), press Run, read
the log, open the report.  Or run from a command line:

    DataChecker.exe <input folder or data workbook> [--out DIR] [--std FILE]
    DataChecker.exe --rerun <..._checked folder>

Output goes to <input folder>_checked/ next to the input:
    REPORT.html, data/ (all XLSX inputs + results), plots/.
"""
from __future__ import annotations

import os
import sys
import threading
import traceback
import webbrowser

# make the bundled tool modules importable both from source and inside PyInstaller
_HERE = getattr(sys, '_MEIPASS', None) or os.path.dirname(os.path.abspath(__file__))
if _HERE not in sys.path:
    sys.path.insert(0, _HERE)

import matplotlib  # noqa: E402
matplotlib.use('Agg')

from lasercal import checker  # noqa: E402


def run_cli(argv):
    import argparse
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('input', nargs='?', help='folder of Qtegra CSVs, or data folder / workbook')
    ap.add_argument('--out', help='output folder (default: <input>_checked)')
    ap.add_argument('--std', help='standards workbook (default: bundled LAICPMS_Standard_Values.xlsx)')
    ap.add_argument('--rerun', help='re-reduce an existing *_checked folder after editing data/*.xlsx')
    ap.add_argument('--no-open', action='store_true', help='do not open the report in a browser')
    a = ap.parse_args(argv)
    if a.rerun:
        report = checker.rerun_folder(a.rerun, a.std)
    elif a.input:
        report = checker.run_checker(a.input, a.out, a.std)
    else:
        ap.print_help(); return 2
    if not a.no_open:
        webbrowser.open('file:///' + report.replace('\\', '/'))
    return 0


def run_gui():
    import tkinter as tk
    from tkinter import filedialog, messagebox, scrolledtext

    root = tk.Tk()
    root.title('LA-ICP-MS Data Checker')
    root.geometry('820x560')
    path_var = tk.StringVar()
    mode_var = tk.StringVar(value='check')

    frm = tk.Frame(root, padx=10, pady=10); frm.pack(fill='both', expand=True)
    tk.Label(frm, text='Input: a folder of Qtegra CSV files, a data folder (workbook + RunOrder.xlsx), '
                       'or an existing *_checked folder to re-run after editing its settings.',
             justify='left', wraplength=780).pack(anchor='w')
    row = tk.Frame(frm); row.pack(fill='x', pady=6)
    tk.Entry(row, textvariable=path_var).pack(side='left', fill='x', expand=True)
    tk.Button(row, text='Browse folder…', command=lambda: path_var.set(filedialog.askdirectory() or path_var.get())).pack(side='left', padx=4)
    tk.Button(row, text='Browse workbook…', command=lambda: path_var.set(
        filedialog.askopenfilename(filetypes=[('Excel', '*.xlsx')]) or path_var.get())).pack(side='left')
    row2 = tk.Frame(frm); row2.pack(fill='x')
    tk.Radiobutton(row2, text='New check (choose parameters automatically)', variable=mode_var, value='check').pack(side='left')
    tk.Radiobutton(row2, text='Re-run a *_checked folder with its edited settings', variable=mode_var, value='rerun').pack(side='left', padx=12)
    log_box = scrolledtext.ScrolledText(frm, height=20, font=('Consolas', 9)); log_box.pack(fill='both', expand=True, pady=6)
    btns = tk.Frame(frm); btns.pack(fill='x')
    run_btn = tk.Button(btns, text='Run', width=12); run_btn.pack(side='left')
    open_btn = tk.Button(btns, text='Open report', width=12, state='disabled'); open_btn.pack(side='left', padx=6)
    state = {'report': None}

    def log(msg):
        log_box.after(0, lambda: (log_box.insert('end', msg + '\n'), log_box.see('end')))

    def work():
        try:
            p = path_var.get().strip().strip('"')
            if not p:
                raise ValueError('Choose an input first.')
            if mode_var.get() == 'rerun':
                rep = checker.rerun_folder(p, log=log)
            else:
                rep = checker.run_checker(p, log=log)
            state['report'] = rep
            log('DONE.')
            root.after(0, lambda: open_btn.config(state='normal'))
        except Exception as ex:
            log('ERROR: ' + str(ex))
            log(traceback.format_exc())
            root.after(0, lambda: messagebox.showerror('Data Checker', str(ex)))
        finally:
            root.after(0, lambda: run_btn.config(state='normal'))

    def on_run():
        run_btn.config(state='disabled'); open_btn.config(state='disabled'); log_box.delete('1.0', 'end')
        threading.Thread(target=work, daemon=True).start()

    run_btn.config(command=on_run)
    open_btn.config(command=lambda: state['report'] and webbrowser.open('file:///' + state['report'].replace('\\', '/')))
    root.mainloop()


if __name__ == '__main__':
    if len(sys.argv) > 1:
        sys.exit(run_cli(sys.argv[1:]))
    run_gui()
