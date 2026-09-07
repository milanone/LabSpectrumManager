import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, RadioButtons, Button
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
import pandas as pd
import os
import tkinter as tk
from tkinterdnd2 import DND_FILES, TkinterDnD

class ScatteringCorrectionApp:
    def __init__(self):
        self.root = TkinterDnD.Tk()
        self.root.title("Scattering Correction - Toolbar Attiva")
        self.root.geometry("1500x900")

        self.x_data = self.y_data = self.y_corr = None
        self.current_file = None

        # --- LAYOUT ---
        # Frame Sinistro: Dati numerici
        self.left_frame = tk.Frame(self.root, width=400, bg="#ffffff")
        self.left_frame.pack(side=tk.LEFT, fill=tk.Y, padx=10, pady=10)
        
        # Frame Destro: Grafico e Toolbar
        self.right_frame = tk.Frame(self.root)
        self.right_frame.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True)

        # Widget Testo con Scrollbar
        tk.Label(self.left_frame, text="nm\tOrig\tCorr", bg="#ffffff", font=('Courier New', 10, 'bold')).pack(anchor='w')
        self.scroll = tk.Scrollbar(self.left_frame)
        self.scroll.pack(side=tk.RIGHT, fill=tk.Y)
        self.data_display = tk.Text(self.left_frame, width=45, font=('Courier New', 9), yscrollcommand=self.scroll.set, wrap='none')
        self.data_display.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        self.scroll.config(command=self.data_display.yview)

        # --- FIGURA MATPLOTLIB ---
        self.fig, self.ax = plt.subplots(figsize=(10, 7))
        plt.subplots_adjust(bottom=0.35) # Spazio per slider e pulsanti
        
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.right_frame)
        self.canvas_widget = self.canvas.get_tk_widget()
        self.canvas_widget.pack(fill=tk.BOTH, expand=True)

        # --- CONTROLLI ZOOM / PAN (TOOLBAR) ---
        # Questa riga aggiunge i tasti per muoversi nello spettro
        self.toolbar = NavigationToolbar2Tk(self.canvas, self.right_frame)
        self.toolbar.update()
        self.toolbar.pack(side=tk.BOTTOM, fill=tk.X)

        # Drag & Drop
        self.canvas_widget.drop_target_register(DND_FILES)
        self.canvas_widget.dnd_bind('<<Drop>>', self.on_drop)

        self._setup_controls()

    def _setup_controls(self):
        # A: Range aggiornato come da tua richiesta (0 - 0.01)
        self.sld_a = Slider(self.fig.add_axes([0.2, 0.22, 0.5, 0.025]), 'A', 0.0, 0.01, valinit=0.0005, valfmt='%.5f')
        
        # P: Esponente (2, 3, 4)
        self.rad_p = RadioButtons(self.fig.add_axes([0.82, 0.12, 0.08, 0.12]), ('2', '3', '4'))
        
        # VS: Esponente logaritmico per la base 10
        self.sld_vs = Slider(self.fig.add_axes([0.2, 0.17, 0.5, 0.025]), 'VS (log)', -2.0, 2.0, valinit=0.0, valfmt='10^%.2f')
        
        # Offset: Traslazione verticale pura
        self.sld_off = Slider(self.fig.add_axes([0.2, 0.12, 0.5, 0.025]), 'Offset', -0.5, 0.5, valinit=0.0, valfmt='%.3f')

        for s in [self.sld_a, self.sld_vs, self.sld_off]:
            s.on_changed(self.update)
        self.rad_p.on_clicked(self.update)

        # Pulsante Esporta
        ax_exp = self.fig.add_axes([0.82, 0.05, 0.1, 0.04])
        self.btn_export = Button(ax_exp, 'Salva CSV', color='#f0f0f0', hovercolor='#ffffff')
        self.btn_export.on_clicked(self.export_data)

    def on_drop(self, event):
        path = event.data.strip('{}')
        if os.path.exists(path):
            try:
                df = pd.read_csv(path)
                self.x_data = pd.to_numeric(df.iloc[:,0], errors='coerce').dropna().values
                self.y_data = pd.to_numeric(df.iloc[:,1], errors='coerce').dropna().values
                self.current_file = path
                
                self.ax.clear()
                self.ax.grid(True, alpha=0.3)
                self.line_orig, = self.ax.plot(self.x_data, self.y_data, color='blue', alpha=0.3, label='Originale')
                self.line_corr, = self.ax.plot(self.x_data, self.x_data*0, color='red', linestyle='--', label='Scattering')
                self.line_res, = self.ax.plot(self.x_data, self.y_data, color='green', lw=1.5, label='Corretto')
                self.ax.legend()
                self.update(None)
            except Exception as e: print(f"Errore: {e}")

    def update(self, val):
        if self.x_data is None: return
        
        A, P, VS, OFF = self.sld_a.val, int(self.rad_p.value_selected), self.sld_vs.val, self.sld_off.val
        
        # FORMULA RICHIESTA: Correzione = (A * 10^VS) * lambda^-P + Offset
        w_norm = self.x_data / 1000.0
        correction = (A * (10**VS)) * (w_norm**-P) + OFF
        
        self.y_corr = self.y_data - correction
        self.line_corr.set_ydata(correction)
        self.line_res.set_ydata(self.y_corr)
        
        # Aggiornamento dati a sinistra (nm, Orig, Corr)
        self.data_display.delete('1.0', tk.END)
        data_lines = [f"{self.x_data[i]:.2f}\t{self.y_data[i]:.4f}\t{self.y_corr[i]:.4f}" for i in range(len(self.x_data))]
        self.data_display.insert(tk.END, "\n".join(data_lines))

        # Autoscale Y
        self.ax.set_ylim(min(np.min(self.y_corr), 0) - 0.05, np.max(self.y_data) + 0.05)
        self.canvas.draw_idle()

    def export_data(self, event):
        if self.y_corr is None: return
        out_path = self.current_file.replace(".csv", "_corrected.csv")
        pd.DataFrame({'nm': self.x_data, 'orig': self.y_data, 'corr': self.y_corr}).to_csv(out_path, index=False)
        self.ax.set_title(f"FILE SALVATO!", color='green', fontsize=10)
        self.canvas.draw_idle()

    def run(self):
        self.root.mainloop()

if __name__ == "__main__":
    app = ScatteringCorrectionApp()
    app.run()