import numpy as np
import matplotlib.pyplot as plt
import os
import sys
import struct
import re
import subprocess
import tkinter as tk
from tkinter import filedialog, scrolledtext, messagebox, Listbox
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
import pandas as pd

try:
    from tkinterdnd2 import DND_FILES, TkinterDnD
    HAS_DND = True
except ImportError:
    HAS_DND = False

CONVERTER_PATH = r"C:\Program Files (x86)\PerkinElmerSP2CSV\PerkinElmerSP2CSV.exe"

class LabSpectrumManager:
    def __init__(self, root):
        self.root = root
        self.root.title("Lab Spectrum Manager (UV-Vis + FTIR + Fluorescence)")
        self.root.geometry("1500x900")

        self.spectra = {}
        self.v_line = None
        self.v_text = None
        self.current_dir = os.getcwd()

        if HAS_DND:
            self.root.drop_target_register(DND_FILES)
            self.root.dnd_bind('<<Drop>>', self.handle_drop)

        # --- INTERFACCIA ---
        self.menu_bar = tk.Menu(root)
        self.file_menu = tk.Menu(self.menu_bar, tearoff=0)
        self.file_menu.add_command(label="Open Files (.dsp, .sp, .spc)", command=self.carica_da_dialog)
        self.file_menu.add_command(label="Export CSV", command=self.esporta_csv)
        self.file_menu.add_separator()
        self.file_menu.add_command(label="Exit", command=root.quit)
        self.menu_bar.add_cascade(label="File", menu=self.file_menu)
        self.root.config(menu=self.menu_bar)

        self.paned = tk.PanedWindow(root, orient=tk.HORIZONTAL, sashrelief=tk.RAISED, sashwidth=4)
        self.paned.pack(fill=tk.BOTH, expand=True)

        # 1. SINISTRA: DATA TABLE
        self.f_data = tk.Frame(self.paned, bg='#ffffff')
        tk.Label(self.f_data, text="DATA TABLE", font=('Arial', 10, 'bold'), bg='#ffffff').pack(pady=5)
        self.data_box = scrolledtext.ScrolledText(self.f_data, width=45, font=('Consolas', 10), bd=0)
        self.data_box.pack(padx=2, pady=2, fill=tk.BOTH, expand=True)
        self.paned.add(self.f_data, width=400)

        # 2. CENTRO: Grafico
        self.f_plot = tk.Frame(self.paned)
        self.fig, self.ax = plt.subplots(figsize=(6, 6))
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.f_plot)
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        self.toolbar = NavigationToolbar2Tk(self.canvas, self.f_plot)
        self.canvas.mpl_connect('motion_notify_event', self.on_mouse_move)
        self.canvas.mpl_connect('axes_leave_event', self.on_mouse_leave)
        self.paned.add(self.f_plot, width=700)

        # 3. DESTRA: Gestione File
        self.f_right = tk.Frame(self.paned, bg='#f0f4f7')
        tk.Label(self.f_right, text="LOADED SPECTRA", font=('Arial', 10, 'bold'), bg='#f0f4f7').pack(pady=5)
        self.file_listbox = Listbox(self.f_right, selectmode=tk.MULTIPLE, height=8, font=('Arial', 9))
        self.file_listbox.pack(padx=10, pady=5, fill=tk.X)

        btn_frame = tk.Frame(self.f_right, bg='#f0f4f7')
        btn_frame.pack(fill=tk.X, padx=10)
        tk.Button(btn_frame, text="Remove Selected", command=self.remove_selected).pack(side=tk.LEFT, expand=True, fill=tk.X, padx=2)
        tk.Button(btn_frame, text="Clear All", command=self.clear_all).pack(side=tk.LEFT, expand=True, fill=tk.X, padx=2)

        tk.Label(self.f_right, text="METADATA", font=('Arial', 10, 'bold'), bg='#f0f4f7').pack(pady=(15, 5))
        self.param_box = scrolledtext.ScrolledText(self.f_right, width=45, font=('Consolas', 9), bg='#f0f4f7', bd=0)
        self.param_box.pack(padx=10, pady=5, fill=tk.BOTH, expand=True)
        self.paned.add(self.f_right, width=400)

        if len(sys.argv) > 1:
            for i in range(1, len(sys.argv)):
                self.processa_file(sys.argv[i])
            self.aggiorna_vista()

    # --- CURSORE ---
    def on_mouse_move(self, event):
        if not event.inaxes or not self.spectra:
            return
        x = event.xdata
        if self.v_line: self.v_line.remove()
        if self.v_text: self.v_text.remove()
        self.v_line = self.ax.axvline(x=x, color='red', linestyle='--', lw=0.7)
        cursor_txt = f"X: {x:.2f}\n"
        for name, data in self.spectra.items():
            df = data['df']
            idx = (df.index.to_series() - x).abs().idxmin()
            cursor_txt += f"{name[:15]}: {df.loc[idx, name]:.4f}\n"
        self.v_text = self.ax.text(
            0.02, 0.98, cursor_txt, transform=self.ax.transAxes,
            verticalalignment='top', fontsize=8,
            bbox=dict(facecolor='white', alpha=0.7)
        )
        self.canvas.draw_idle()

    def on_mouse_leave(self, event):
        if self.v_line: self.v_line.remove(); self.v_line = None
        if self.v_text: self.v_text.remove(); self.v_text = None
        self.canvas.draw_idle()

    # --- FILE I/O ---
    def handle_drop(self, event):
        files = self.root.tk.splitlist(event.data)
        for f in files:
            self.processa_file(f.strip('{}'))
        self.current_dir = os.path.dirname(files[-1].strip('{}'))
        self.aggiorna_vista()

    def carica_da_dialog(self):
        paths = filedialog.askopenfilenames(initialdir=self.current_dir,
                                            filetypes=[("Spectra Files", "*.dsp *.sp *.csv *.spc *.SPC"), ("All Files", "*.*")])
        if paths:
            for p in paths: self.processa_file(p)
            self.current_dir = os.path.dirname(paths[-1])
            self.aggiorna_vista()

    def processa_file(self, path):
        ext = os.path.splitext(path)[1].lower()
        try:
            if ext == '.dsp':
                self.leggi_dsp(path)
            elif ext == '.sp':
                self.leggi_sp(path)
            elif ext == '.csv':
                self.leggi_csv(path)
            elif ext == '.spc':
                self.leggi_spc(path)
        except Exception as e:
            messagebox.showerror("Error", f"Could not read {os.path.basename(path)}: {e}")

    def leggi_dsp(self, path):
        with open(path, 'r', encoding='latin-1') as f:
            lines = [line.strip() for line in f.readlines()]
        name = os.path.basename(path).replace('.dsp', '')
        wl_start, wl_end, n_punti = float(lines[5]), float(lines[6]), int(lines[8])
        idx_data = lines.index("#DATA") + 1
        y_vals = [float(lines[i]) for i in range(idx_data, idx_data + n_punti)]
        x_vals = np.linspace(wl_start, wl_end, n_punti)
        self.spectra[name] = {
            'df': pd.DataFrame({'x': x_vals, name: y_vals}).set_index('x'),
            'info': {'Sample': name, 'Date': lines[16], 'Type': 'UV-Vis'}
        }

    def leggi_sp(self, path):
        subprocess.run([CONVERTER_PATH, path], check=True, creationflags=0x08000000)
        csv_path = path.replace('.sp', '.sp.csv').replace('.SP', '.SP.csv')
        if os.path.exists(csv_path):
            df_temp = pd.read_csv(csv_path)
            name = os.path.basename(path).replace('.sp', '')
            x_vals = df_temp.iloc[:, 0].values
            y_vals = df_temp.iloc[:, 1].values
            self.spectra[name] = {
                'df': pd.DataFrame({'x': x_vals, name: y_vals}).set_index('x'),
                'info': {'Sample': name, 'Date': 'N/A', 'Type': 'FTIR'}
            }
            # CSV temporaneo mantenuto intenzionalmente (utile per debug)

    def leggi_csv(self, path):
        df = pd.read_csv(path, index_col=0)
        tipo = 'FTIR' if df.index.max() > 2000 else 'UV-Vis'
        for col in df.columns:
            serie = df[[col]].dropna()
            self.spectra[col] = {
                'df': serie,
                'info': {'Sample': col, 'Date': 'N/A', 'Type': tipo}
            }

    def leggi_spc(self, path):
        """Parser per file binari .SPC (Shimadzu RF-5301 fluorescence spectrometer)."""
        with open(path, 'rb') as f:
            raw = f.read()

        name = os.path.splitext(os.path.basename(path))[0]
        header_text = raw[:512].decode('latin-1', errors='ignore')

        # Dati numerici
        x_start, x_end = struct.unpack('<ff', raw[79:87])
        y_segment = raw[235:]
        n_points = len(y_segment) // 4
        y_data = np.frombuffer(y_segment[:n_points * 4], dtype=np.float32).copy()
        x_data = np.linspace(x_start, x_end, len(y_data))
        pitch = abs(x_data[1] - x_data[0]) if len(x_data) > 1 else 0.0

        # Wavelength fissa (eccitazione o emissione)
        wl_fissa = 0.0
        for i in range(40, 160):
            try:
                val = struct.unpack('<f', raw[i:i + 4])[0]
                if 200.0 <= val <= 900.0 and val % 1 == 0:
                    wl_fissa = val
                    break
            except:
                continue

        # Fenditure (Ex e Em slit width)
        slits = []
        for i in range(40, 180):
            try:
                val = struct.unpack('<f', raw[i:i + 4])[0]
                if val in [1.5, 3.0, 5.0, 10.0, 15.0, 20.0]:
                    slits.append(val)
            except:
                continue

        # Metadati testuali
        modo = "EM" if "EM" in header_text else "EX" if "EX" in header_text else "N/D"
        campione = raw[154:210].decode('latin-1', errors='ignore').split('\x00')[0].strip()
        if not campione:
            campione = name

        # Data/ora
        area_data = raw[210:245].decode('latin-1', errors='ignore')
        ora = re.search(r'(\d{2}:\d{2})', area_data)
        giorno = re.search(r'(\d{2}/\d{2}/\d{2})', area_data)
        data_str = f"{ora.group(1) if ora else 'N/D'} {giorno.group(1) if giorno else 'N/D'}"

        scan_range = f"{min(x_data):.1f} - {max(x_data):.1f} nm"
        wl_label = "Excitation WL" if modo == "EM" else "Emission WL"

        self.spectra[name] = {
            'df': pd.DataFrame({'x': x_data, name: y_data}).set_index('x'),
            'info': {
                'Sample':          campione,
                'Date':            data_str,
                'Type':            'Fluorescence',
                'Instrument':      'Shimadzu RF-5301',
                'Scan Type':       modo,
                'Scan Range':      scan_range,
                'Sample Pitch':    f"{pitch:.1f} nm",
                'Ex Slit Width':   f"{slits[0] if len(slits) > 0 else 5.0} nm",
                'Em Slit Width':   f"{slits[1] if len(slits) > 1 else 5.0} nm",
                wl_label:          f"{wl_fissa:.1f} nm",
                'Scan Speed':      next((s for s in ["Very Fast", "Fast", "Medium", "Slow", "Super"] if s in header_text), "N/D"),
                'Sensitivity':     "High" if "High" in header_text else "Low",
                'Response Time':   "Auto" if "Auto" in header_text else "N/D",
                'Shutter':         ("Manual" if "Manual" in header_text else "Auto") +
                                   " / " + ("Open" if "Open" in header_text else "Closed"),
            }
        }

    # --- VISTA ---
    def aggiorna_vista(self):
        self.file_listbox.delete(0, tk.END)
        self.data_box.delete(1.0, tk.END)
        self.param_box.delete(1.0, tk.END)
        self.ax.clear()
        self.v_line = self.v_text = None

        if not self.spectra:
            self.canvas.draw()
            return

        for sname in self.spectra.keys():
            self.file_listbox.insert(tk.END, sname)

        master_df = pd.concat([s['df'] for s in self.spectra.values()], axis=1).sort_index()

        # Tabella tabulata
        self.data_box.insert(tk.END, master_df.to_csv(sep='\t', lineterminator='\n'))

        # Plotting — rilevamento tipo dal campo 'Type'
        types = {d['info']['Type'] for d in self.spectra.values()}

        for col in master_df.columns:
            valid = master_df[col].dropna()
            self.ax.plot(valid.index, valid.values, label=col, lw=1.5)

        self.ax.set_xlim(master_df.index.min(), master_df.index.max())

        if 'FTIR' in types:
            self.ax.set_xlabel("Wavenumber (cm⁻¹)")
            self.ax.set_ylabel("Transmittance / Absorbance")
            self.ax.invert_xaxis()
        elif 'Fluorescence' in types:
            self.ax.set_xlabel("Wavelength (nm)")
            self.ax.set_ylabel("Intensity (A.U.)")
        else:
            self.ax.set_xlabel("Wavelength (nm)")
            self.ax.set_ylabel("Absorbance (A)")

        self.ax.legend(fontsize='8', loc='best')
        self.ax.grid(True, linestyle=':', alpha=0.6)
        self.canvas.draw()

        # Metadata panel — display generico di tutte le chiavi info
        for name, data in self.spectra.items():
            self.param_box.insert(tk.END, f"{'='*38}\n")
            for k, v in data['info'].items():
                self.param_box.insert(tk.END, f"{k:<20}: {v}\n")
            self.param_box.insert(tk.END, "\n")

    # --- GESTIONE LISTA ---
    def remove_selected(self):
        selected = self.file_listbox.curselection()
        for i in reversed(selected):
            name = self.file_listbox.get(i)
            if name in self.spectra: del self.spectra[name]
        self.aggiorna_vista()

    def clear_all(self):
        self.spectra = {}
        self.aggiorna_vista()

    def esporta_csv(self):
        if not self.spectra: return
        path = filedialog.asksaveasfilename(defaultextension=".csv", initialdir=self.current_dir)
        if path:
            dfs = [s['df'] for s in self.spectra.values()]
            pd.concat(dfs, axis=1).to_csv(path)

if __name__ == "__main__":
    root = TkinterDnD.Tk() if HAS_DND else tk.Tk()
    app = LabSpectrumManager(root)
    root.mainloop()
