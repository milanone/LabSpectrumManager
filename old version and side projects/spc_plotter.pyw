import numpy as np
import matplotlib.pyplot as plt
import struct
import re
import csv
import os
import tkinter as tk
from tkinter import filedialog, scrolledtext, messagebox
from tkinterdnd2 import DND_FILES, TkinterDnD
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk

class ShimadzuAnalyzer:
    def __init__(self, root):
        self.root = root
        self.root.title("Shimadzu SPC Expert Pro")
        self.root.geometry("1200x800")
        
        # Stato dell'applicazione
        self.x_data = None
        self.y_data = None
        self.current_path = ""
        self.current_filename = ""
        self.current_dir = os.getcwd() # Cartella di avvio predefinita

        # Configurazione Drag and Drop
        self.root.drop_target_register(DND_FILES)
        self.root.dnd_bind('<<Drop>>', self.handle_drop)

        # Barra dei Menu
        self.menu_bar = tk.Menu(root)
        self.file_menu = tk.Menu(self.menu_bar, tearoff=0)
        self.file_menu.add_command(label="Open SPC", command=self.carica_da_dialog)
        self.file_menu.add_command(label="Export to CSV", command=self.esporta_csv)
        self.file_menu.add_separator()
        self.file_menu.add_command(label="Exit", command=root.quit)
        self.menu_bar.add_cascade(label="File", menu=self.file_menu)
        self.root.config(menu=self.menu_bar)

        # Layout Principale
        self.f_left = tk.Frame(root)
        self.f_left.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        
        self.f_right = tk.Frame(root, width=450, bg='#fdf5e6')
        self.f_right.pack(side=tk.RIGHT, fill=tk.Y, padx=10, pady=10)
        
        tk.Label(self.f_right, text="SPECTRUM PARAMETERS", font=('Arial', 10, 'bold'), bg='#fdf5e6').pack(pady=5)
        
        self.box = scrolledtext.ScrolledText(self.f_right, width=55, height=35, font=('Consolas', 10))
        self.box.pack(padx=5, pady=5, fill=tk.BOTH, expand=True)
        self.box.insert(tk.END, ">>> Trascina qui un file .SPC oppure usa File -> Open")

        # Setup Grafico
        self.fig, self.ax = plt.subplots(figsize=(8, 6))
        self.ax.set_facecolor('#fcfcfc')
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.f_left)
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        self.toolbar = NavigationToolbar2Tk(self.canvas, self.f_left)
        self.toolbar.update()

    def handle_drop(self, event):
        path = event.data.strip('{}') # Rimuove le graffe per percorsi con spazi
        if path.lower().endswith('.spc'):
            self.processa_file(path)
        else:
            messagebox.showwarning("Formato non valido", "Trascina solo file .SPC")

    def carica_da_dialog(self):
        path = filedialog.askopenfilename(
            initialdir=self.current_dir,
            title="Seleziona file Shimadzu SPC",
            filetypes=[("File SPC", "*.SPC")]
        )
        if path:
            self.processa_file(path)

    def processa_file(self, path):
        try:
            # Sincronizzazione Directory
            self.current_path = os.path.abspath(path)
            self.current_dir = os.path.dirname(self.current_path)
            self.current_filename = os.path.basename(self.current_path)
            os.chdir(self.current_dir)

            with open(self.current_path, 'rb') as f:
                raw_data = f.read()

            header_text = raw_data[:512].decode('latin-1', errors='ignore')

            # --- ESTRAZIONE DATI NUMERICI ---
            x_start, x_end = struct.unpack('<ff', raw_data[79:87])
            y_segment = raw_data[235:]
            n_points = len(y_segment) // 4
            self.y_data = np.frombuffer(y_segment[:n_points*4], dtype=np.float32).copy()
            self.x_data = np.linspace(x_start, x_end, len(self.y_data))
            
            # Calcoli derivati
            pitch = abs(self.x_data[1] - self.x_data[0]) if len(self.x_data) > 1 else 0.0
            scan_range = f"{min(self.x_data):.1f} - {max(self.x_data):.1f} nm"

            # --- AGGIORNAMENTO GRAFICO ---
            self.ax.clear()
            self.ax.plot(self.x_data, self.y_data, color='#2ecc71', lw=2.5) # Verde Calceina
            self.ax.set_title(f"Spectrum: {self.current_filename}", fontweight='bold', pad=15)
            self.ax.set_xlabel("Wavelength (nm)")
            self.ax.set_ylabel("Intensity (A.U.)")
            self.ax.grid(True, linestyle=':', alpha=0.5)
            self.canvas.draw()

            # --- ESTRAZIONE PARAMETRI (LOGICA VALIDATA) ---
            wl_fissa = 0.0
            for i in range(40, 160):
                try:
                    val = struct.unpack('<f', raw_data[i:i+4])[0]
                    if 200.0 <= val <= 900.0 and val % 1 == 0:
                        wl_fissa = val
                        break
                except: continue

            slits = []
            for i in range(40, 180):
                try:
                    val = struct.unpack('<f', raw_data[i:i+4])[0]
                    if val in [1.5, 3.0, 5.0, 10.0, 15.0, 20.0]: slits.append(val)
                except: continue

            modo = "EM" if "EM" in header_text else "EX" if "EX" in header_text else "N/D"
            campione = raw_data[154:210].decode('latin-1', errors='ignore').split('\x00')[0].strip()
            
            # Costruzione Report
            res =  f"CAMPIONE: {campione}\n" + "="*50 + "\n"
            params = [
                ("Strumento", "Shimadzu RF-5301"),
                ("Soft. Ver.", "1.1"),
                ("Scan Type", modo),
                ("Scan Speed", next((s for s in ["Very Fast", "Fast", "Medium", "Slow", "Super"] if s in header_text), "N/D")),
                ("Scan Range", scan_range),
                ("Sensitivity", "High" if "High" in header_text else "Low"),
                ("Response Time", "Auto" if "Auto" in header_text else "N/D"),
                ("Ex Slit Width", f"{slits[0] if len(slits)>0 else 5.0} nm"),
                ("Em Slit Width", f"{slits[1] if len(slits)>1 else 5.0} nm"),
                ("Sample Pitch", f"{pitch:.1f} nm"),
                ("Shutter Mode", "Manual" if "Manual" in header_text else "Auto"),
                ("Shutter Position", "Open" if "Open" in header_text else "Closed"),
                ("Date/Time", self.estrai_data_string(raw_data)),
                ("Excitation WL" if modo == "EM" else "Emission WL", f"{wl_fissa:.1f} nm")
            ]
            for k, v in params:
                res += f"{k:<24}: {v}\n"
            
            res += "="*50 + f"\nPunti totali: {len(self.y_data)}"
            
            self.box.delete(1.0, tk.END)
            self.box.insert(tk.END, res)

        except Exception as e:
            messagebox.showerror("Errore", f"Impossibile analizzare il file:\n{str(e)}")

    def estrai_data_string(self, dati):
        area = dati[210:245].decode('latin-1', errors='ignore')
        o = re.search(r'(\d{2}:\d{2})', area)
        d = re.search(r'(\d{2}/\d{2}/\d{2})', area)
        return f"{o.group(1) if o else 'N/D'} {d.group(1) if d else 'N/D'}"

    def esporta_csv(self):
        if self.x_data is None:
            messagebox.showwarning("Export", "Nessun dato caricato da esportare.")
            return
        
        path_csv = filedialog.asksaveasfilename(
            initialdir=self.current_dir,
            defaultextension=".csv", 
            initialfile=f"{self.current_filename[:-4]}_data.csv",
            title="Esporta dati spettro in CSV",
            filetypes=[("CSV files", "*.csv"), ("All files", "*.*")]
        )
        
        if path_csv:
            try:
                # Se l'utente salva in una cartella diversa, la memorizziamo per la prossima volta
                self.current_dir = os.path.dirname(path_csv)
                
                with open(path_csv, 'w', newline='') as f:
                    writer = csv.writer(f)
                    writer.writerow(['Wavelength (nm)', 'Intensity (A.U.)'])
                    for row in zip(self.x_data, self.y_data):
                        writer.writerow([f"{row[0]:.3f}", f"{row[1]:.4f}"])
                
                messagebox.showinfo("Export", "File CSV generato con successo!")
            except Exception as e:
                messagebox.showerror("Errore Export", f"Errore durante il salvataggio:\n{e}")

if __name__ == "__main__":
    # Avvio con supporto Drag & Drop
    root = TkinterDnD.Tk()
    app = ShimadzuAnalyzer(root)
    root.mainloop()