import numpy as np
import matplotlib.pyplot as plt
import os
import sys
import struct
import re
import datetime
import pickle
import traceback
import tkinter as tk
from tkinter import filedialog, scrolledtext, messagebox, Listbox
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
import pandas as pd

try:
    from tkinterdnd2 import DND_FILES, TkinterDnD
    HAS_DND = True
except ImportError:
    HAS_DND = False

# Stile "Origin-like" condiviso col repo fratello PlotStyleKit (così LabSpectrumManager,
# KleistekManager e pca-gui usano lo stesso modulo senza duplicarlo). Caricato per path,
# cercando PRIMA una copia locale in questa cartella (per fissare una versione specifica
# per questo solo programma, se mai servisse) e POI la cartella fratella ../PlotStyleKit
# condivisa; fallback no-op se manca del tutto (vedi il warning mostrato all'avvio in
# __init__). I rcParams vengono impostati subito, prima di creare qualunque figura, così
# i grafici nascono già in stile Origin.
def _carica_origin_style():
    import importlib.util as ilu
    here = os.path.dirname(os.path.abspath(__file__))
    for path in (os.path.join(here, 'origin_style.py'),
                 os.path.join(here, '..', 'PlotStyleKit', 'origin_style.py')):
        if os.path.isfile(path):
            spec = ilu.spec_from_file_location('origin_style', path)
            mod = ilu.module_from_spec(spec)
            spec.loader.exec_module(mod)
            return mod
    return None

try:
    origin_style = _carica_origin_style()
    if origin_style is not None:
        origin_style.applica_rcparams()
except Exception:
    origin_style = None


class LabSpectrumManager:
    def __init__(self, root):
        self.root = root
        self.root.title("Lab Spectrum Manager (UV-Vis + FTIR + Fluorescence)")
        # Dimensione iniziale commisurata allo schermo (non un fisso 1500x900, che su
        # schermi piccoli o con scaling elevato sporge sotto la taskbar): al più il
        # 90% della larghezza e l'80% dell'altezza disponibili (l'80%, non l'85%, per
        # lasciare margine alla barra delle applicazioni), finestra centrata.
        screen_w = self.root.winfo_screenwidth()
        screen_h = self.root.winfo_screenheight()
        win_w = min(1500, int(screen_w * 0.90))
        win_h = min(900, int(screen_h * 0.80))
        pos_x = (screen_w - win_w) // 2
        pos_y = (screen_h - win_h) // 2
        self.root.geometry(f"{win_w}x{win_h}+{pos_x}+{pos_y}")

        if origin_style is None:
            # avviso non bloccante: l'app funziona comunque (stile matplotlib di
            # default), ma senza il look Origin né "Edit Figure..."; lo si nota
            # solo guardando il grafico se non si avvisa esplicitamente qui.
            self.root.after(200, lambda: messagebox.showwarning(
                "PlotStyleKit non trovato",
                "Il repo condiviso PlotStyleKit (origin_style.py / plot_editor.pyw) non è stato "
                "trovato né come copia locale in questa cartella né come cartella fratella "
                "..\\PlotStyleKit.\n\n"
                "I grafici useranno lo stile matplotlib di default (niente look Origin) e "
                "\"Edit Figure...\" non sarà disponibile.\n\n"
                "Clona https://github.com/milanone/PlotStyleKit accanto a questo progetto per "
                "abilitarli."))

        self._pe_module = None       # modulo plot_editor, caricato alla prima necessità

        self.spectra = {}
        self._hidden_spectra = set()   # nomi nascosti dal grafico (restano in lista/tabella)
        self._dirty = False   # True se ci sono risultati calcolati non ancora esportati
        self.v_line = None
        self.v_text = None
        self.current_dir = os.getcwd()
        self._cursor_data = []   # cache (name, x_array, y_array) per il cursore
        self._cursor_box_pos = None   # (ax_x, ax_y, ha, va) se trascinato manualmente; None = posizione "best" automatica
        self._cursor_drag = False     # True mentre l'utente trascina il riquadro X/Y

        # Stato pannello scattering
        self._sc_name     = None
        self._sc_x        = None
        self._sc_y        = None
        self._sc_lines    = {}
        self._sc_vline_lo = None
        self._sc_vline_hi = None
        self._sc_drag     = None   # 'lo' | 'hi' | None
        self._sc_cids     = []     # canvas connection IDs
        self._sc_lim_cids   = []   # callback id per xlim/ylim_changed
        self._sc_user_zoomed = False   # True se l'utente ha zoomato/pannato
        self._sc_setting_lim = False   # guard: stiamo impostando noi i limiti

        # Stato pannello scaled subtraction
        self._sub_name_a = None
        self._sub_name_b = None
        self._sub_index_a = None  # pandas Index di A (usato per il DataFrame risultato)
        self._sub_x      = None   # asse x di A come numpy array (per il grafico)
        self._sub_y_a    = None
        self._sub_y_b    = None
        self._sub_lines  = {}

        # Stato pannello normalize
        self._norm_names = []
        self._norm_lines = {}

        # Stato pannello baseline FTIR
        self._bl_name    = None
        self._bl_x       = None
        self._bl_y       = None
        self._bl_lines   = {}
        self._bl_va      = None   # linea ancora 1 (trascinabile)
        self._bl_vb      = None   # linea ancora 2 (trascinabile)
        self._bl_drag    = None   # 'a' | 'b' | None
        self._bl_cids    = []
        self._bl_syncing = False

        # Stato pannello baseline adattiva
        self._ab_name  = None
        self._ab_x     = None
        self._ab_y     = None
        self._ab_lines = {}

        # Stato pannello smoothing
        self._sm_name  = None
        self._sm_x     = None
        self._sm_y     = None
        self._sm_lines = {}

        # Stato pannello deconvoluzione (lorentziane)
        self._dc_name   = None
        self._dc_x      = None
        self._dc_y      = None
        self._dc_peaks  = []      # lista di dict {'x0','A','w'}
        self._dc_offset = 0.0
        self._dc_fitted = False
        self._dc_cids   = []
        self._dc_range  = None    # [lo, hi] intervallo di fit
        self._dc_va     = None    # linea intervallo sinistra
        self._dc_vb     = None    # linea intervallo destra
        self._dc_drag   = None    # 'a' | 'b' | None
        self._dc_drag_moved = False   # distingue un click fermo (→ aggiungi picco) da un drag reale

        # Stato pannello trim (ritaglio di una porzione di spettro)
        self._tr_name     = None
        self._tr_x        = None
        self._tr_y        = None
        self._tr_lines    = {}
        self._tr_vline_lo = None
        self._tr_vline_hi = None
        self._tr_drag     = None   # 'lo' | 'hi' | None
        self._tr_cids     = []     # canvas connection IDs

        if HAS_DND:
            self.root.drop_target_register(DND_FILES)
            self.root.dnd_bind('<<Drop>>', self.handle_drop)

        # --- MENU ---
        self.menu_bar = tk.Menu(root)
        self.file_menu = tk.Menu(self.menu_bar, tearoff=0)
        self.file_menu.add_command(label="Open Files (.dsp, .sp, .spc, .csv)", command=self.carica_da_dialog)
        self.file_menu.add_command(label="Export CSV", command=self.esporta_csv)
        self.file_menu.add_command(label="Save Figure (pickle)", command=self.salva_figura_pickle)
        self.file_menu.add_command(label="Edit Figure...", command=self.apri_editor_figura)
        self.file_menu.add_separator()
        self.file_menu.add_command(label="Exit", command=self.on_exit)
        self.menu_bar.add_cascade(label="File", menu=self.file_menu)
        self.root.config(menu=self.menu_bar)
        self.root.protocol("WM_DELETE_WINDOW", self.on_exit)

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
        # Il margine superiore di default di matplotlib (12% dell'altezza) è pensato
        # per un titolo del grafico che qui non è mai impostato: ridotto per non
        # lasciare una fascia vuota sopra il grafico. Impostato una sola volta sulla
        # figura condivisa: resta valido per tutti i pannelli (ax.clear() non lo tocca).
        self.fig.subplots_adjust(top=0.96)
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.f_plot)
        self.toolbar = NavigationToolbar2Tk(self.canvas, self.f_plot, pack_toolbar=False)
        self.toolbar.pack(side=tk.TOP, fill=tk.X)
        self.canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)
        self.canvas.mpl_connect('motion_notify_event', self.on_mouse_move)
        self.canvas.mpl_connect('axes_leave_event', self.on_mouse_leave)
        self.canvas.mpl_connect('button_press_event', self._cursor_drag_start)
        self.canvas.mpl_connect('button_release_event', self._cursor_drag_stop)
        # Incolla dati x,y dalla clipboard: clic sul grafico per il focus, poi Ctrl+V
        plot_widget = self.canvas.get_tk_widget()
        plot_widget.bind('<Button-1>', lambda e: plot_widget.focus_set(), add='+')
        plot_widget.bind('<Control-v>', self._incolla_clipboard)
        plot_widget.bind('<Control-V>', self._incolla_clipboard)
        # add='+' è essenziale: senza, questa bind sostituirebbe quella interna di
        # FigureCanvasTk su <Button-3> (che genera il button_press_event di matplotlib),
        # e il tasto destro smetterebbe di funzionare in deconvoluzione (e ovunque)
        # tranne che con un doppio click, l'unico ancora instradato correttamente
        plot_widget.bind('<Button-3>', self._menu_grafico, add='+')
        self.paned.add(self.f_plot, width=700)

        # 3. DESTRA: Gestione File
        self.f_right = tk.Frame(self.paned, bg='#f0f4f7')
        tk.Label(self.f_right, text="LOADED SPECTRA", font=('Arial', 10, 'bold'), bg='#f0f4f7').pack(pady=5)
        self.file_listbox = Listbox(self.f_right, selectmode=tk.MULTIPLE, height=8, font=('Arial', 9))
        self.file_listbox.pack(padx=10, pady=5, fill=tk.X)
        self.file_listbox.bind('<Button-3>', self._listbox_tasto_destro)

        btn_frame = tk.Frame(self.f_right, bg='#f0f4f7')
        btn_frame.pack(fill=tk.X, padx=10)
        tk.Button(btn_frame, text="Remove Selected", command=self.remove_selected).pack(side=tk.LEFT, expand=True, fill=tk.X, padx=2)
        tk.Button(btn_frame, text="Clear All", command=self.clear_all).pack(side=tk.LEFT, expand=True, fill=tk.X, padx=2)

        btn_frame_vis = tk.Frame(self.f_right, bg='#f0f4f7')
        btn_frame_vis.pack(fill=tk.X, padx=10, pady=(2, 0))
        tk.Button(btn_frame_vis, text="Hide Selected", command=self.hide_selected).pack(side=tk.LEFT, expand=True, fill=tk.X, padx=2)
        tk.Button(btn_frame_vis, text="Show Selected", command=self.show_selected).pack(side=tk.LEFT, expand=True, fill=tk.X, padx=2)

        btn_frame2 = tk.Frame(self.f_right, bg='#f0f4f7')
        btn_frame2.pack(fill=tk.X, padx=10, pady=(2, 0))
        tk.Button(btn_frame2, text="Average Selected", command=self.media_selezionati,
                  bg='#e8f5e9', font=('Arial', 9)).pack(side=tk.LEFT, expand=True, fill=tk.X, padx=2)
        tk.Button(btn_frame2, text="Scaled Subtraction", command=self.apri_sottrazione_scalata,
                  bg='#fff3e0', font=('Arial', 9)).pack(side=tk.LEFT, expand=True, fill=tk.X, padx=2)

        tk.Button(self.f_right, text="Scattering Correction",
                  command=self.apri_correzione_scattering,
                  bg='#dce8f5', font=('Arial', 9)).pack(fill=tk.X, padx=10, pady=(4, 0))
        tk.Button(self.f_right, text="Normalize Selected",
                  command=self.apri_normalizzazione,
                  bg='#f3e5f5', font=('Arial', 9)).pack(fill=tk.X, padx=10, pady=(4, 0))
        tk.Button(self.f_right, text="Baseline (FTIR)",
                  command=self.apri_baseline,
                  bg='#d7ccc8', font=('Arial', 9)).pack(fill=tk.X, padx=10, pady=(4, 0))
        tk.Button(self.f_right, text="Adaptive Baseline",
                  command=self.apri_baseline_adattiva,
                  bg='#c8e6c9', font=('Arial', 9)).pack(fill=tk.X, padx=10, pady=(4, 0))
        tk.Button(self.f_right, text="Smooth Selected",
                  command=self.apri_smoothing,
                  bg='#ffe0b2', font=('Arial', 9)).pack(fill=tk.X, padx=10, pady=(4, 0))
        tk.Button(self.f_right, text="Deconvolution",
                  command=self.apri_deconvoluzione,
                  bg='#b3e5fc', font=('Arial', 9)).pack(fill=tk.X, padx=10, pady=(4, 0))
        tk.Button(self.f_right, text="Trim",
                  command=self.apri_trim,
                  bg='#cfd8dc', font=('Arial', 9)).pack(fill=tk.X, padx=10, pady=(4, 0))

        # --- PANNELLO METADATA (visibile per default) ---
        self.f_metadata = tk.Frame(self.f_right, bg='#f0f4f7')
        self.f_metadata.pack(fill=tk.BOTH, expand=True)
        tk.Label(self.f_metadata, text="METADATA", font=('Arial', 10, 'bold'), bg='#f0f4f7').pack(pady=(15, 5))
        self.param_box = scrolledtext.ScrolledText(self.f_metadata, width=45, font=('Consolas', 9), bg='#f0f4f7', bd=0)
        self.param_box.pack(padx=10, pady=5, fill=tk.BOTH, expand=True)

        # --- PANNELLO SCATTERING (nascosto per default) ---
        self.f_scattering = tk.Frame(self.f_right, bg='#f0f4f7')
        # non packed — appare solo quando attivato
        self._costruisci_pannello_scattering()

        # --- PANNELLO SCALED SUBTRACTION (nascosto per default) ---
        self.f_subtraction = tk.Frame(self.f_right, bg='#f0f4f7')
        # non packed — appare solo quando attivato
        self._costruisci_pannello_sottrazione()

        # --- PANNELLO NORMALIZE (nascosto per default) ---
        self.f_normalize = tk.Frame(self.f_right, bg='#f0f4f7')
        # non packed — appare solo quando attivato
        self._costruisci_pannello_normalizzazione()

        # --- PANNELLO BASELINE FTIR (nascosto per default) ---
        self.f_baseline = tk.Frame(self.f_right, bg='#f0f4f7')
        # non packed — appare solo quando attivato
        self._costruisci_pannello_baseline()

        # --- PANNELLO BASELINE ADATTIVA (nascosto per default) ---
        self.f_adaptive = tk.Frame(self.f_right, bg='#f0f4f7')
        # non packed — appare solo quando attivato
        self._costruisci_pannello_adattivo()

        # --- PANNELLO SMOOTHING (nascosto per default) ---
        self.f_smooth = tk.Frame(self.f_right, bg='#f0f4f7')
        # non packed — appare solo quando attivato
        self._costruisci_pannello_smoothing()

        # --- PANNELLO DECONVOLUZIONE (nascosto per default) ---
        self.f_deconv = tk.Frame(self.f_right, bg='#f0f4f7')
        # non packed — appare solo quando attivato
        self._costruisci_pannello_deconv()

        # --- PANNELLO TRIM (nascosto per default) ---
        self.f_trim = tk.Frame(self.f_right, bg='#f0f4f7')
        # non packed — appare solo quando attivato
        self._costruisci_pannello_trim()

        self.paned.add(self.f_right, width=400)

        if len(sys.argv) > 1:
            for i in range(1, len(sys.argv)):
                self.processa_file(sys.argv[i])
            self.aggiorna_vista()

    # --- PANNELLO SCATTERING: costruzione widget (una sola volta) ---
    def _costruisci_pannello_scattering(self):
        bg = '#f0f4f7'

        tk.Label(self.f_scattering, text="SCATTERING CORRECTION",
                 font=('Arial', 10, 'bold'), bg=bg).pack(pady=(15, 2))
        self._sc_label_nome = tk.Label(self.f_scattering, text="", font=('Consolas', 9, 'italic'), bg=bg)
        self._sc_label_nome.pack()

        ctrl = tk.Frame(self.f_scattering, bg=bg)
        ctrl.pack(fill=tk.X, padx=10, pady=6)

        # Tre parametri, ciascuno slider fine + casella editabile (valori a piacimento)
        self._sc_VS  = self._make_linked_row(ctrl, 'VS',        0.0,    0.3,   0.0002,   '{:.4f}')
        self._sc_M   = self._make_linked_row(ctrl, 'Slope (m)', -0.001, 0.001, 0.000002, '{:.6f}')
        self._sc_M_lock = tk.BooleanVar(value=True)
        m_lock_row = tk.Frame(ctrl, bg=bg)
        m_lock_row.pack(fill=tk.X, pady=(0, 3))
        tk.Label(m_lock_row, text='', width=9, bg=bg).pack(side=tk.LEFT)   # allineamento
        tk.Checkbutton(m_lock_row, text='fix slope', variable=self._sc_M_lock, bg=bg,
                       font=('Consolas', 8)).pack(side=tk.LEFT, padx=(10, 0))
        self._sc_OFF = self._make_linked_row(ctrl, 'Offset',    -1.5,   1.5,   0.001,    '{:.3f}')

        # Esponente P: slider continuo + casella, con preset 1/2/3/4 e blocco "fix P"
        self._sc_P = self._make_linked_row(ctrl, 'P', 1.0, 5.0, 0.01, '{:.2f}')
        self._linked_set(self._sc_P, 2.0)
        preset = tk.Frame(ctrl, bg=bg)
        preset.pack(fill=tk.X, pady=(0, 3))
        tk.Label(preset, text='', width=9, bg=bg).pack(side=tk.LEFT)   # allineamento
        for p_val in (1, 2, 3, 4):
            tk.Button(preset, text=str(p_val), width=2, font=('Consolas', 8),
                      command=lambda v=p_val: self._linked_set(self._sc_P, float(v))).pack(side=tk.LEFT, padx=2)
        self._sc_P_lock = tk.BooleanVar(value=False)
        tk.Checkbutton(preset, text='fix P', variable=self._sc_P_lock, bg=bg,
                       font=('Consolas', 8)).pack(side=tk.LEFT, padx=(10, 0))

        # Range label + pulsante Auto
        self._sc_range_label = tk.Label(self.f_scattering, text="Fit range: —",
                                        font=('Consolas', 8), bg='#f0f4f7', fg='#555555')
        self._sc_range_label.pack(pady=(6, 0))
        tk.Button(self.f_scattering, text="Auto (fit window)",
                  command=self._sc_auto, bg='#fff9c4', font=('Arial', 9)).pack(fill=tk.X, padx=10, pady=(2, 0))

        # Pulsanti Apply / Cancel
        btn_sc = tk.Frame(self.f_scattering, bg=bg)
        btn_sc.pack(fill=tk.X, padx=10, pady=10)
        tk.Button(btn_sc, text="Apply", command=self._sc_applica,
                  bg='#2ecc71', fg='white', font=('Arial', 9, 'bold')).pack(side=tk.LEFT, expand=True, fill=tk.X, padx=2)
        tk.Button(btn_sc, text="Cancel", command=self._sc_annulla).pack(side=tk.LEFT, expand=True, fill=tk.X, padx=2)

    # --- SCATTERING: riga parametro con slider fine + casella editabile ---
    def _make_linked_row(self, parent, label, lo, hi, res, fmt):
        """Crea una riga slider+Entry sincronizzati. La casella accetta valori
        a piacimento (anche oltre il range dello slider); lo slider resta clampato.
        Restituisce un dict di stato con chiavi 'value','slider','entry','lo','hi','fmt'."""
        bg = '#f0f4f7'
        st = {'syncing': False, 'lo': lo, 'hi': hi, 'fmt': fmt,
              'value':  tk.DoubleVar(value=0.0),
              'slider': tk.DoubleVar(value=0.0),
              'entry':  tk.StringVar(value=fmt.format(0.0))}
        row = tk.Frame(parent, bg=bg)
        row.pack(fill=tk.X, pady=3)
        tk.Label(row, text=label, width=9, anchor='w', bg=bg, font=('Consolas', 9)).pack(side=tk.LEFT)
        sc = tk.Scale(row, variable=st['slider'], from_=lo, to=hi, resolution=res,
                      orient=tk.HORIZONTAL, length=185, bg=bg, font=('Consolas', 8),
                      showvalue=False,
                      command=lambda v, s=st: self._linked_from_slider(s, v))
        sc.pack(side=tk.LEFT)
        st['scale'] = sc
        e = tk.Entry(row, textvariable=st['entry'], width=8, font=('Consolas', 9))
        e.pack(side=tk.LEFT, padx=4)
        e.bind('<Return>',   lambda ev, s=st: self._linked_from_entry(s))
        e.bind('<FocusOut>', lambda ev, s=st: self._linked_from_entry(s))
        st['value'].trace_add('write', self._sc_aggiorna_preview)
        return st

    def _linked_set_range(self, st, lo, hi, res):
        """Ridimensiona il range dello slider (la casella resta comunque libera)."""
        st['lo'], st['hi'] = lo, hi
        st['scale'].config(from_=lo, to=hi, resolution=res)

    def _linked_set(self, st, value, origine=None):
        """Imposta il parametro allineando valore reale, slider (clampato) e casella."""
        st['syncing'] = True
        if origine != 'entry':
            st['entry'].set(st['fmt'].format(value))
        if origine != 'slider':
            st['slider'].set(min(max(value, st['lo']), st['hi']))
        st['value'].set(value)   # trace → preview
        st['syncing'] = False

    def _linked_from_slider(self, st, val):
        if st['syncing']:
            return
        self._linked_set(st, float(val), origine='slider')

    def _linked_from_entry(self, st):
        if st['syncing']:
            return
        try:
            v = float(st['entry'].get().replace(',', '.'))
        except ValueError:
            return
        self._linked_set(st, v, origine='entry')

    # --- SCATTERING: apertura pannello ---
    def apri_correzione_scattering(self):
        sel = self._selezione_effettiva()
        if len(sel) != 1:
            messagebox.showwarning("Scattering Correction", "Seleziona esattamente uno spettro dalla lista.")
            return
        name = self.file_listbox.get(sel[0])
        if self.spectra[name]['info']['Type'] == 'FTIR':
            messagebox.showwarning("Scattering Correction", "La correzione scattering non è applicabile a spettri FTIR.")
            return

        self._sc_name = name
        self._sc_x = self.spectra[name]['df'].index.to_numpy(dtype=float)
        self._sc_y = self.spectra[name]['df'][name].to_numpy(dtype=float)

        # Range dei parametri commisurati all'ampiezza Y (e X per la pendenza):
        # offset e VS ~ ampiezza dello spettro; slope ~ ampiezza/estensione X
        span = float(self._sc_y.max() - self._sc_y.min())
        if span <= 0:
            span = max(abs(float(self._sc_y.max())), 0.1)
        x_absmax = max(abs(float(self._sc_x.min())), abs(float(self._sc_x.max())), 1.0)
        m_hi = span / x_absmax
        # VS è il valore a 1000; la correzione all'estremo di bassa λ è
        # VS·(λmin/1000)^-2. Scaliamo VS così che tale correzione arrivi a ~2·ampiezza,
        # tenendo conto dell'amplificazione (altrimenti il range risulta troppo ampio).
        lam_min = float(self._sc_x.min())
        vs_hi = 2.0 * span * (lam_min / 1000.0) ** 2 if lam_min > 0 else span
        self._linked_set_range(self._sc_VS,  0.0,    vs_hi, vs_hi / 1000.0)
        self._linked_set_range(self._sc_M,  -m_hi,   m_hi,  2 * m_hi / 1000.0)
        self._linked_set_range(self._sc_OFF, -span,  span,  2 * span / 1000.0)

        # Reset slider ai valori di default
        self._linked_set(self._sc_VS,  0.0)
        self._linked_set(self._sc_M,   0.0)
        self._sc_M_lock.set(True)
        self._linked_set(self._sc_OFF, 0.0)
        self._linked_set(self._sc_P, 2.0)
        self._sc_label_nome.config(text=name)

        # Prepara grafico principale con le tre linee di anteprima
        self.ax.clear()
        self.v_line = self.v_text = None
        tipo = self.spectra[name]['info']['Type']
        self.ax.set_xlabel("Raman shift (cm⁻¹)" if tipo == 'Raman' else "Wavelength (nm)")
        self.ax.set_ylabel("Intensity (A.U.)" if tipo in ('Fluorescence', 'Raman') else "Absorbance (A)")
        self.ax.grid(True, linestyle=':', alpha=0.6)

        self._sc_lines['orig'], = self.ax.plot(self._sc_x, self._sc_y,
                                               color='steelblue', alpha=0.5, lw=1.5, label=f'{name} (orig)')
        self._sc_lines['corr'], = self.ax.plot(self._sc_x, np.zeros_like(self._sc_y),
                                               color='orange', lw=1.2, ls='--', label='Scattering')
        self._sc_lines['res'],  = self.ax.plot(self._sc_x, self._sc_y.copy(),
                                               color='green', lw=2, label=f'{name}_corr')
        self.ax.legend(fontsize=8)

        # Linee verticali trascinabili per la finestra di fit: di default coprono
        # l'intero spettro (l'auto-fit vincolato è ormai robusto su tutta la banda)
        x_lo = self._sc_x.min()
        x_hi = self._sc_x.max()
        self._sc_vline_lo = self.ax.axvline(x=x_lo, color='purple', ls='--', lw=1.5, alpha=0.8)
        self._sc_vline_hi = self.ax.axvline(x=x_hi, color='purple', ls='--', lw=1.5, alpha=0.8)
        self._sc_drag = None
        self._sc_cids = [
            self.canvas.mpl_connect('button_press_event',   self._sc_drag_start),
            self.canvas.mpl_connect('motion_notify_event',  self._sc_drag_move),
            self.canvas.mpl_connect('button_release_event', self._sc_drag_stop),
        ]
        self._sc_aggiorna_range_label()
        self.canvas.draw()

        # Rileva zoom/pan dell'utente per non sovrascriverlo durante la preview.
        # Connesso dopo il primo draw così l'autoscale iniziale non conta come zoom.
        self._sc_user_zoomed = False
        self._sc_lim_cids = [
            self.ax.callbacks.connect('xlim_changed', self._sc_on_lim_changed),
            self.ax.callbacks.connect('ylim_changed', self._sc_on_lim_changed),
        ]

        # Mostra pannello scattering, nascondi metadata
        self.f_metadata.pack_forget()
        self.f_scattering.pack(fill=tk.BOTH, expand=True)

    def _sc_on_lim_changed(self, ax=None):
        # Ignora i cambi di limiti che impostiamo noi (guard); gli altri = zoom utente
        if not self._sc_setting_lim:
            self._sc_user_zoomed = True

    # --- SCATTERING: aggiornamento preview live ---
    def _sc_aggiorna_preview(self, *_):
        if self._sc_x is None or 'corr' not in self._sc_lines:
            return
        VS, P, M, OFF = (self._sc_VS['value'].get(), self._sc_P['value'].get(),
                         self._sc_M['value'].get(), self._sc_OFF['value'].get())
        correction = VS * (self._sc_x / 1000.0) ** -P - M * self._sc_x + OFF
        y_corr = self._sc_y - correction
        self._sc_lines['corr'].set_ydata(correction)
        self._sc_lines['res'].set_ydata(y_corr)
        # Autoscale solo finché l'utente non ha zoomato: così lo zoom viene preservato
        if not self._sc_user_zoomed:
            y_all  = np.concatenate([self._sc_y, y_corr, correction])
            margin = (y_all.max() - y_all.min()) * 0.05 or 0.05
            self._sc_setting_lim = True
            self.ax.set_ylim(y_all.min() - margin, y_all.max() + margin)
            self._sc_setting_lim = False
        self.canvas.draw_idle()

    # --- SCATTERING: fit automatico sulla finestra definita dalle linee ---
    def _sc_auto(self):
        if self._sc_x is None or self._sc_vline_lo is None:
            return
        P0 = float(self._sc_P['value'].get())
        lock_P = self._sc_P_lock.get()
        M0_fixed = float(self._sc_M['value'].get())
        lock_M = self._sc_M_lock.get()
        x, y = self._sc_x, self._sc_y
        x_lo = self._sc_vline_lo.get_xdata()[0]
        x_hi = self._sc_vline_hi.get_xdata()[0]
        x_lo, x_hi = min(x_lo, x_hi), max(x_lo, x_hi)
        mask = (x >= x_lo) & (x <= x_hi)
        if mask.sum() < 3:
            messagebox.showwarning("Auto fit", "La finestra selezionata contiene meno di 3 punti.")
            return
        x_ref, y_ref = x[mask], y[mask]

        # Fit lineare iniziale a P fisso: VS*(x/1000)^-P - M*x + OFF = y (sulla finestra)
        A_ref = np.column_stack([
            (x_ref / 1000.0) ** -P0,
            -x_ref,
            np.ones_like(x_ref),
        ])
        if lock_M:
            coeffs, _, _, _ = np.linalg.lstsq(A_ref[:, [0, 2]], y_ref + M0_fixed * x_ref, rcond=None)
            VS0, M0, OFF0 = float(max(coeffs[0], 0.0)), M0_fixed, float(coeffs[1])
        else:
            coeffs, _, _, _ = np.linalg.lstsq(A_ref, y_ref, rcond=None)
            VS0, M0, OFF0 = float(max(coeffs[0], 0.0)), float(coeffs[1]), float(coeffs[2])

        # Punto di partenza ammissibile: abbassa l'offset finché la traccia è sotto lo spettro
        corr_full = VS0 * (x / 1000.0) ** -P0 - M0 * x + OFF0
        over = float(np.max(corr_full - y))
        if over > 0:
            OFF0 -= over

        VS_fit, M_fit, OFF_fit, P_fit = VS0, M0, OFF0, P0
        step = max(1, len(x) // 800)          # decima i vincoli (la traccia è liscia)
        xc, yc = x[::step], y[::step]
        # Scale per condizionare l'ottimizzazione (VS, M, x su ordini diversi)
        s = np.maximum(np.abs(A_ref).max(axis=0), 1e-12)

        # Rifinitura vincolata (traccia ≤ spettro in ogni punto). Richiede scipy;
        # senza scipy resta la soluzione ammissibile con solo l'offset abbassato.
        # Se lo slope è fissato, il suo bound collassa a un punto (M0*s[1]) e SLSQP
        # lo tiene bloccato pur restando nella stessa formulazione del fit libero
        m_bound = (M0 * s[1], M0 * s[1]) if lock_M else (None, None)

        try:
            from scipy.optimize import minimize
            if lock_P:
                # P fisso → modello lineare nei parametri: usa il jacobiano analitico
                Ac = np.column_stack([(xc / 1000.0) ** -P0, -xc, np.ones_like(xc)])
                Ars, Acs = A_ref / s, Ac / s
                cons = [{'type': 'ineq',
                         'fun': lambda q: yc - Acs @ q,
                         'jac': lambda q: -Acs}]
                res = minimize(lambda q: float((Ars @ q - y_ref) @ (Ars @ q - y_ref)),
                               np.array([VS0, M0, OFF0]) * s,
                               jac=lambda q: 2.0 * Ars.T @ (Ars @ q - y_ref),
                               method='SLSQP',
                               bounds=[(0.0, None), m_bound, (None, None)],
                               constraints=cons, options={'maxiter': 500, 'ftol': 1e-12})
                if np.all(np.isfinite(res.x)):
                    VS_fit, M_fit, OFF_fit = (res.x / s)
            else:
                # P libero → modello non lineare: ottimizza anche l'esponente in [1, 5]
                def corr(q, xx):
                    vs, m, off, pp = q[0] / s[0], q[1] / s[1], q[2] / s[2], q[3]
                    return vs * (xx / 1000.0) ** -pp - m * xx + off
                cons = [{'type': 'ineq', 'fun': lambda q: yc - corr(q, xc)}]
                q0 = np.array([VS0 * s[0], M0 * s[1], OFF0 * s[2], P0])
                res = minimize(lambda q: float((corr(q, x_ref) - y_ref) @ (corr(q, x_ref) - y_ref)),
                               q0, method='SLSQP',
                               bounds=[(0.0, None), m_bound, (None, None), (1.0, 5.0)],
                               constraints=cons, options={'maxiter': 800, 'ftol': 1e-12})
                if np.all(np.isfinite(res.x)):
                    VS_fit = res.x[0] / s[0]
                    M_fit  = res.x[1] / s[1]
                    OFF_fit = res.x[2] / s[2]
                    P_fit  = float(res.x[3])
        except ImportError:
            pass

        VS_fit = float(max(VS_fit, 0.0))

        # Garanzia finale su TUTTI i punti (il vincolo SLSQP era su griglia decimata)
        corr_full = VS_fit * (x / 1000.0) ** -P_fit - M_fit * x + OFF_fit
        over = float(np.max(corr_full - y))
        if over > 0:
            OFF_fit -= over

        # Usa i valori reali del fit (VS ≥ 0); se eccedono il range dello slider
        # la casella li mostra comunque e lo slider resta agganciato al massimo
        self._linked_set(self._sc_VS,  round(VS_fit, 6))
        self._linked_set(self._sc_OFF, round(float(OFF_fit), 6))
        if not lock_M:
            self._linked_set(self._sc_M, round(float(M_fit), 8))
        if not lock_P:
            self._linked_set(self._sc_P, round(P_fit, 3))

    # --- SCATTERING: drag delle linee verticali ---
    def _sc_aggiorna_range_label(self):
        if self._sc_vline_lo is None:
            return
        lo = self._sc_vline_lo.get_xdata()[0]
        hi = self._sc_vline_hi.get_xdata()[0]
        lo, hi = min(lo, hi), max(lo, hi)
        self._sc_range_label.config(text=f"Fit range: {lo:.1f} – {hi:.1f}")

    def _sc_drag_start(self, event):
        if event.inaxes != self.ax or event.button != 1:
            return
        for key, line in [('lo', self._sc_vline_lo), ('hi', self._sc_vline_hi)]:
            x_disp = self.ax.transData.transform((line.get_xdata()[0], 0))[0]
            if abs(event.x - x_disp) < 8:
                self._sc_drag = key
                return

    def _sc_drag_move(self, event):
        if self._sc_drag is None or event.inaxes != self.ax or event.xdata is None:
            return
        xv = event.xdata
        lo = self._sc_vline_lo.get_xdata()[0]
        hi = self._sc_vline_hi.get_xdata()[0]
        if self._sc_drag == 'lo':
            xv = min(xv, hi - 1)
            self._sc_vline_lo.set_xdata([xv, xv])
        else:
            xv = max(xv, lo + 1)
            self._sc_vline_hi.set_xdata([xv, xv])
        self._sc_aggiorna_range_label()
        self.canvas.draw_idle()

    def _sc_drag_stop(self, event):
        self._sc_drag = None

    # --- SCATTERING: applica e torna alla vista normale ---
    def _sc_applica(self):
        VS, P, M, OFF = (self._sc_VS['value'].get(), self._sc_P['value'].get(),
                         self._sc_M['value'].get(), self._sc_OFF['value'].get())
        correction = VS * (self._sc_x / 1000.0) ** -P - M * self._sc_x + OFF
        y_corr   = self._sc_y - correction
        new_name = f"{self._sc_name}_corr"
        new_info = dict(self.spectra[self._sc_name]['info'])
        new_info['Sample']         = new_name
        new_info['Derived from']   = self._sc_name
        new_info['Scatter VS']     = f"{VS:.4f}"
        new_info['Scatter P']      = str(P)
        new_info['Scatter Slope']  = f"{M:.6f}"
        new_info['Scatter Offset'] = f"{OFF:.3f}"
        self.spectra[new_name] = {
            'df':   pd.DataFrame({'x': self._sc_x, new_name: y_corr}).set_index('x'),
            'info': new_info,
        }
        self._dirty = True
        self._sc_chiudi_pannello()
        self.aggiorna_vista()

    # --- SCATTERING: annulla senza modifiche ---
    def _sc_annulla(self):
        self._sc_chiudi_pannello()
        self.aggiorna_vista()

    # --- SCATTERING: chiude il pannello e ripristina metadata ---
    def _sc_chiudi_pannello(self):
        for cid in self._sc_cids:
            self.canvas.mpl_disconnect(cid)
        self._sc_cids     = []
        for cid in self._sc_lim_cids:
            self.ax.callbacks.disconnect(cid)
        self._sc_lim_cids   = []
        self._sc_user_zoomed = False
        self._sc_drag     = None
        self._sc_vline_lo = None
        self._sc_vline_hi = None
        self._sc_name     = None
        self._sc_x        = None
        self._sc_y        = None
        self._sc_lines    = {}
        self.f_scattering.pack_forget()
        self.f_metadata.pack(fill=tk.BOTH, expand=True)

    # --- PANNELLO SCALED SUBTRACTION: costruzione widget ---
    def _costruisci_pannello_sottrazione(self):
        bg = '#f0f4f7'

        tk.Label(self.f_subtraction, text="SCALED SUBTRACTION",
                 font=('Arial', 10, 'bold'), bg=bg).pack(pady=(15, 2))

        # Etichetta con i nomi degli spettri coinvolti
        self._sub_label_ab = tk.Label(self.f_subtraction, text="", font=('Consolas', 9, 'italic'),
                                      bg=bg, wraplength=360, justify='left')
        self._sub_label_ab.pack(padx=10)

        ctrl = tk.Frame(self.f_subtraction, bg=bg)
        ctrl.pack(fill=tk.X, padx=10, pady=8)

        row = tk.Frame(ctrl, bg=bg)
        row.pack(fill=tk.X, pady=4)
        tk.Label(row, text='k', width=5, anchor='w', bg=bg, font=('Consolas', 10, 'bold')).pack(side=tk.LEFT)
        self._var_k = tk.DoubleVar(value=1.0)
        tk.Scale(row, variable=self._var_k, from_=-5.0, to=5.0, resolution=0.001,
                 orient=tk.HORIZONTAL, length=270, bg=bg,
                 font=('Consolas', 8), showvalue=True).pack(side=tk.LEFT)

        # Spettro di riferimento
        ref_row = tk.Frame(ctrl, bg=bg)
        ref_row.pack(fill=tk.X, pady=4)
        tk.Label(ref_row, text='Ref', width=5, anchor='w', bg=bg, font=('Consolas', 10)).pack(side=tk.LEFT)
        self._sub_ref_var = tk.StringVar(value='— none —')
        self._sub_ref_menu = tk.OptionMenu(ref_row, self._sub_ref_var, '— none —')
        self._sub_ref_menu.config(font=('Consolas', 8), width=28)
        self._sub_ref_menu.pack(side=tk.LEFT, fill=tk.X, expand=True)
        self._sub_ref_var.trace_add('write', self._sub_aggiorna_riferimento)

        # Formula display
        self._sub_formula = tk.Label(self.f_subtraction, text="Result = A − 1.000 × B",
                                     font=('Consolas', 9), bg=bg, fg='#555555')
        self._sub_formula.pack()

        self._var_k.trace_add('write', self._sub_aggiorna_preview)

        # Pulsante Swap
        tk.Button(self.f_subtraction, text="Swap A ↔ B", command=self._sub_scambia_ab,
                  bg='#e8eaf6', font=('Arial', 9)).pack(fill=tk.X, padx=10, pady=(0, 4))

        # Pulsanti Apply / Cancel
        btn_sub = tk.Frame(self.f_subtraction, bg=bg)
        btn_sub.pack(fill=tk.X, padx=10, pady=10)
        tk.Button(btn_sub, text="Apply", command=self._sub_applica,
                  bg='#2ecc71', fg='white', font=('Arial', 9, 'bold')).pack(side=tk.LEFT, expand=True, fill=tk.X, padx=2)
        tk.Button(btn_sub, text="Cancel", command=self._sub_annulla).pack(side=tk.LEFT, expand=True, fill=tk.X, padx=2)

    # --- SCALED SUBTRACTION: apertura pannello ---
    def apri_sottrazione_scalata(self):
        sel = self.file_listbox.curselection()
        if len(sel) != 2:
            messagebox.showwarning("Scaled Subtraction", "Seleziona esattamente 2 spettri (A e B, in quest'ordine).")
            return
        name_a = self.file_listbox.get(sel[0])
        name_b = self.file_listbox.get(sel[1])

        types = {self.spectra[name_a]['info']['Type'], self.spectra[name_b]['info']['Type']}
        if len(types) > 1:
            messagebox.showwarning("Scaled Subtraction", f"Gli spettri hanno tipi diversi: {', '.join(types)}.")
            return

        self._sub_name_a = name_a
        self._sub_name_b = name_b

        # Asse x di A; B allineato per indice (stessa x garantita dall'utente)
        df_a = self.spectra[name_a]['df']
        df_b = self.spectra[name_b]['df']
        self._sub_index_a = df_a.index
        self._sub_x   = df_a.index.to_numpy(dtype=float)
        self._sub_y_a = df_a[name_a].to_numpy(dtype=float)
        self._sub_y_b = df_b[name_b].reindex(df_a.index).to_numpy(dtype=float)

        self._var_k.set(1.0)
        self._sub_label_ab.config(text=f"A: {name_a}\nB: {name_b}")

        # Aggiorna opzioni menu riferimento
        opts = ['— none —'] + [n for n in self.spectra if n not in (name_a, name_b)]
        m = self._sub_ref_menu['menu']
        m.delete(0, 'end')
        for opt in opts:
            m.add_command(label=opt, command=lambda v=opt: self._sub_ref_var.set(v))
        self._sub_ref_var.set('— none —')

        # Prepara grafico
        self.ax.clear()
        self.v_line = self.v_text = None
        tipo = self.spectra[name_a]['info']['Type']
        if tipo == 'FTIR':
            self.ax.set_xlabel("Wavenumber (cm⁻¹)")
            self.ax.set_ylabel(self._ftir_ylabel([self._sub_name_a, self._sub_name_b]))
        elif tipo == 'Raman':
            self.ax.set_xlabel("Raman shift (cm⁻¹)")
            self.ax.set_ylabel("Intensity (A.U.)")
        elif tipo == 'Fluorescence':
            self.ax.set_xlabel("Wavelength (nm)")
            self.ax.set_ylabel("Intensity (A.U.)")
        else:
            self.ax.set_xlabel("Wavelength (nm)")
            self.ax.set_ylabel("Absorbance (A)")
        self.ax.grid(True, linestyle=':', alpha=0.6)

        self._sub_lines['ref'],  = self.ax.plot([], [], color='gray', lw=1.2,
                                                ls=':', visible=False, label='_nolegend_')
        self._sub_lines['a'],    = self.ax.plot(self._sub_x, self._sub_y_a,
                                                color='steelblue', alpha=0.6, lw=1.5, label=f'{name_a} (A)')
        self._sub_lines['b'],    = self.ax.plot(self._sub_x, self._sub_y_b,
                                                color='tomato', alpha=0.3, lw=1.2, ls='--', label=f'{name_b} (B)')
        self._sub_lines['kb'],   = self.ax.plot(self._sub_x, self._sub_y_b,
                                                color='tomato', alpha=0.8, lw=1.5, label='k×B')
        self._sub_lines['res'],  = self.ax.plot(self._sub_x, self._sub_y_a - self._sub_y_b,
                                                color='green', lw=2, label='A − k×B')
        self.ax.legend(fontsize=8)
        if tipo == 'FTIR':
            self.ax.invert_xaxis()
        self.canvas.draw()

        self.f_metadata.pack_forget()
        self.f_subtraction.pack(fill=tk.BOTH, expand=True)
        self.root.bind('<Left>',  lambda e: self._sub_freccia(-0.001))
        self.root.bind('<Right>', lambda e: self._sub_freccia(+0.001))

    # --- SCALED SUBTRACTION: preview live ---
    def _sub_aggiorna_preview(self, *_):
        if self._sub_x is None:
            return
        k = self._var_k.get()
        kb = k * self._sub_y_b
        result = self._sub_y_a - kb
        self._sub_lines['kb'].set_ydata(kb)
        self._sub_lines['res'].set_ydata(result)
        self._sub_formula.config(text=f"Result = A − {k:.3f} × B")
        y_all  = np.concatenate([self._sub_y_a, self._sub_y_b, kb, result])
        y_min, y_max = np.nanmin(y_all), np.nanmax(y_all)
        margin = (y_max - y_min) * 0.05 or 0.05
        self.ax.set_ylim(y_min - margin, y_max + margin)
        self.canvas.draw_idle()

    # --- SCALED SUBTRACTION: applica ---
    def _sub_applica(self):
        k      = self._var_k.get()
        result = self._round_sig(self._sub_y_a - k * self._sub_y_b, sig=6)
        label  = f"{self._sub_name_a}_minus_{k:.3f}x{self._sub_name_b}"
        # Evita collisioni
        final_label = label
        counter = 2
        while final_label in self.spectra:
            final_label = f"{label}_{counter}"
            counter += 1
        tipo = self.spectra[self._sub_name_a]['info']['Type']
        result_df = pd.DataFrame({final_label: result}, index=self._sub_index_a)
        self.spectra[final_label] = {
            'df': result_df,
            'info': {
                'Sample':    final_label,
                'Date':      'N/A',
                'Type':      tipo,
                'Operation': 'Scaled Subtraction',
                'A':         self._sub_name_a,
                'B':         self._sub_name_b,
                'k':         f"{k:.4f}",
                'Formula':   f"A − {k:.4f} × B",
            }
        }
        self._dirty = True
        self._sub_chiudi_pannello()
        self.aggiorna_vista()

    # --- SCALED SUBTRACTION: spettro di riferimento ---
    def _sub_aggiorna_riferimento(self, *_):
        line = self._sub_lines.get('ref')
        if line is None:
            return
        name = self._sub_ref_var.get()
        if name == '— none —' or name not in self.spectra:
            line.set_visible(False)
            line.set_label('_nolegend_')
        else:
            df = self.spectra[name]['df']
            line.set_xdata(df.index.to_numpy(dtype=float))
            line.set_ydata(df[name].to_numpy(dtype=float))
            line.set_visible(True)
            line.set_label(f'{name} (ref)')
        self.ax.legend(fontsize=8)
        self.canvas.draw_idle()

    # --- SCALED SUBTRACTION: aggiustamento fine con frecce ---
    def _sub_freccia(self, delta):
        if self._sub_x is None:
            return
        nuovo = round(self._var_k.get() + delta, 3)
        self._var_k.set(max(-5.0, min(5.0, nuovo)))

    # --- SCALED SUBTRACTION: inverti A e B ---
    def _sub_scambia_ab(self):
        if self._sub_x is None:
            return
        self._sub_name_a, self._sub_name_b = self._sub_name_b, self._sub_name_a
        df_a = self.spectra[self._sub_name_a]['df']
        df_b = self.spectra[self._sub_name_b]['df']
        self._sub_index_a = df_a.index
        self._sub_x   = df_a.index.to_numpy(dtype=float)
        self._sub_y_a = df_a[self._sub_name_a].to_numpy(dtype=float)
        self._sub_y_b = df_b[self._sub_name_b].reindex(df_a.index).to_numpy(dtype=float)
        self._sub_label_ab.config(text=f"A: {self._sub_name_a}\nB: {self._sub_name_b}")
        self._sub_lines['a'].set_xdata(self._sub_x)
        self._sub_lines['a'].set_ydata(self._sub_y_a)
        self._sub_lines['a'].set_label(f'{self._sub_name_a} (A)')
        self._sub_lines['b'].set_xdata(self._sub_x)
        self._sub_lines['b'].set_ydata(self._sub_y_b)
        self._sub_lines['b'].set_label(f'{self._sub_name_b} (B)')
        self._sub_lines['kb'].set_xdata(self._sub_x)
        self._sub_lines['res'].set_xdata(self._sub_x)
        self.ax.legend(fontsize=8)
        self._sub_aggiorna_preview()

    # --- SCALED SUBTRACTION: annulla ---
    def _sub_annulla(self):
        self._sub_chiudi_pannello()
        self.aggiorna_vista()

    # --- SCALED SUBTRACTION: chiude il pannello ---
    def _sub_chiudi_pannello(self):
        self.root.unbind('<Left>')
        self.root.unbind('<Right>')
        self._sub_name_a  = None
        self._sub_name_b  = None
        self._sub_index_a = None
        self._sub_x       = None
        self._sub_y_a     = None
        self._sub_y_b     = None
        self._sub_lines   = {}
        self.f_subtraction.pack_forget()
        self.f_metadata.pack(fill=tk.BOTH, expand=True)

    # --- PANNELLO NORMALIZE: costruzione widget ---
    def _costruisci_pannello_normalizzazione(self):
        bg = '#f0f4f7'

        tk.Label(self.f_normalize, text="NORMALIZE",
                 font=('Arial', 10, 'bold'), bg=bg).pack(pady=(15, 2))

        self._norm_label_nomi = tk.Label(self.f_normalize, text="",
                                         font=('Consolas', 9, 'italic'),
                                         bg=bg, wraplength=360, justify='left')
        self._norm_label_nomi.pack(padx=10)

        ctrl = tk.Frame(self.f_normalize, bg=bg)
        ctrl.pack(fill=tk.X, padx=10, pady=8)

        row = tk.Frame(ctrl, bg=bg)
        row.pack(fill=tk.X, pady=4)
        tk.Label(row, text='λ norm', width=7, anchor='w', bg=bg,
                 font=('Consolas', 10)).pack(side=tk.LEFT)
        self._norm_wl_var = tk.StringVar()
        tk.Entry(row, textvariable=self._norm_wl_var, width=10,
                 font=('Consolas', 10)).pack(side=tk.LEFT, padx=4)
        self._norm_unit_label = tk.Label(row, text='nm', font=('Consolas', 9),
                                         bg=bg, fg='#555555')
        self._norm_unit_label.pack(side=tk.LEFT)

        self._norm_wl_var.trace_add('write', self._norm_aggiorna_preview)

        self._norm_info_label = tk.Label(self.f_normalize, text="",
                                         font=('Consolas', 8), bg=bg, fg='#555555',
                                         wraplength=360, justify='left')
        self._norm_info_label.pack(padx=10, pady=(0, 4))

        btn_norm = tk.Frame(self.f_normalize, bg=bg)
        btn_norm.pack(fill=tk.X, padx=10, pady=10)
        tk.Button(btn_norm, text="Apply", command=self._norm_applica,
                  bg='#2ecc71', fg='white', font=('Arial', 9, 'bold')).pack(side=tk.LEFT, expand=True, fill=tk.X, padx=2)
        tk.Button(btn_norm, text="Cancel", command=self._norm_annulla).pack(side=tk.LEFT, expand=True, fill=tk.X, padx=2)

    # --- NORMALIZE: apertura pannello ---
    def apri_normalizzazione(self):
        sel = self._selezione_effettiva()
        if not sel:
            messagebox.showwarning("Normalize", "Seleziona almeno uno spettro.")
            return
        self._norm_names = [self.file_listbox.get(i) for i in sel]
        types = {self.spectra[n]['info']['Type'] for n in self._norm_names}
        if len(types) > 1:
            messagebox.showwarning("Normalize", f"Gli spettri selezionati hanno tipi diversi: {', '.join(types)}.")
            return

        tipo = next(iter(types))
        self._norm_label_nomi.config(text=', '.join(self._norm_names))
        self._norm_unit_label.config(text='cm⁻¹' if tipo == 'FTIR' else 'nm')

        # Prepara grafico
        self.ax.clear()
        self.v_line = self.v_text = None
        if tipo == 'FTIR':
            self.ax.set_xlabel("Wavenumber (cm⁻¹)")
            self.ax.set_ylabel(self._ftir_ylabel(self._norm_names))
            self.ax.invert_xaxis()
        elif tipo == 'Raman':
            self.ax.set_xlabel("Raman shift (cm⁻¹)")
            self.ax.set_ylabel("Intensity (A.U.)")
        elif tipo == 'Fluorescence':
            self.ax.set_xlabel("Wavelength (nm)")
            self.ax.set_ylabel("Intensity (A.U.)")
        else:
            self.ax.set_xlabel("Wavelength (nm)")
            self.ax.set_ylabel("Normalized Absorbance")
        self.ax.grid(True, linestyle=':', alpha=0.6)

        colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
        for i, name in enumerate(self._norm_names):
            color = colors[i % len(colors)]
            df = self.spectra[name]['df']
            x = df.index.to_numpy(dtype=float)
            y = df[name].to_numpy(dtype=float)
            self._norm_lines[f'{name}_orig'], = self.ax.plot(
                x, y, color=color, alpha=0.25, lw=1.2, ls='--')
            self._norm_lines[f'{name}_norm'], = self.ax.plot(
                x, y, color=color, lw=1.8, label=f'{name}_norm')

        default_wl = round((self.spectra[self._norm_names[0]]['df'].index.min() +
                            self.spectra[self._norm_names[0]]['df'].index.max()) / 2, 1)
        self._norm_lines['vline'] = self.ax.axvline(
            x=default_wl, color='red', ls=':', lw=1.0, alpha=0.7)
        self.ax.legend(fontsize=8)
        self.canvas.draw()

        self.f_metadata.pack_forget()
        self.f_normalize.pack(fill=tk.BOTH, expand=True)
        self._norm_wl_var.set(str(default_wl))

    # --- NORMALIZE: preview live ---
    def _norm_aggiorna_preview(self, *_):
        if not self._norm_names or 'vline' not in self._norm_lines:
            return
        try:
            wl = float(self._norm_wl_var.get())
        except ValueError:
            return

        info_lines = []
        all_y = []
        for name in self._norm_names:
            df = self.spectra[name]['df']
            x = df.index.to_numpy(dtype=float)
            y = df[name].to_numpy(dtype=float)
            idx = np.argmin(np.abs(x - wl))
            val = y[idx]
            if val == 0:
                info_lines.append(f"{name[:18]}: valore=0, saltato")
                continue
            y_norm = y / val
            self._norm_lines[f'{name}_norm'].set_ydata(y_norm)
            all_y.extend(y_norm[np.isfinite(y_norm)])
            info_lines.append(f"{name[:18]}: {val:.4f} @ {x[idx]:.1f}")

        self._norm_lines['vline'].set_xdata([wl, wl])
        self._norm_info_label.config(text='\n'.join(info_lines))

        if all_y:
            mn, mx = min(all_y), max(all_y)
            margin = (mx - mn) * 0.05 or 0.05
            self.ax.set_ylim(mn - margin, mx + margin)
        self.canvas.draw_idle()

    # --- NORMALIZE: applica e salva ---
    def _norm_applica(self):
        try:
            wl = float(self._norm_wl_var.get())
        except ValueError:
            messagebox.showwarning("Normalize", "Inserisci un valore numerico per la lunghezza d'onda.")
            return
        tipo = self.spectra[self._norm_names[0]]['info']['Type']
        for name in self._norm_names:
            df = self.spectra[name]['df']
            x = df.index.to_numpy(dtype=float)
            y = df[name].to_numpy(dtype=float)
            idx = np.argmin(np.abs(x - wl))
            val = y[idx]
            if val == 0:
                messagebox.showwarning("Normalize",
                    f"'{name}': valore 0 a {x[idx]:.1f}, spettro saltato.")
                continue
            y_norm = self._round_sig(y / val, sig=6)
            new_name = f"{name}_norm{wl:.0f}"
            final_name = new_name
            counter = 2
            while final_name in self.spectra:
                final_name = f"{new_name}_{counter}"
                counter += 1
            self.spectra[final_name] = {
                'df': pd.DataFrame({final_name: y_norm}, index=df.index),
                'info': {
                    'Sample':     final_name,
                    'Date':       'N/A',
                    'Type':       tipo,
                    'Operation':  'Normalize',
                    'Source':     name,
                    'Norm at':    f"{x[idx]:.2f}",
                    'Norm value': f"{val:.6f}",
                }
            }
            self._dirty = True
        self._norm_chiudi_pannello()
        self.aggiorna_vista()

    # --- NORMALIZE: annulla ---
    def _norm_annulla(self):
        self._norm_chiudi_pannello()
        self.aggiorna_vista()

    # --- NORMALIZE: chiude il pannello ---
    def _norm_chiudi_pannello(self):
        self._norm_names = []
        self._norm_lines = {}
        self.f_normalize.pack_forget()
        self.f_metadata.pack(fill=tk.BOTH, expand=True)

    # --- PANNELLO BASELINE FTIR: costruzione widget ---
    def _costruisci_pannello_baseline(self):
        bg = '#f0f4f7'

        tk.Label(self.f_baseline, text="BASELINE (FTIR)",
                 font=('Arial', 10, 'bold'), bg=bg).pack(pady=(15, 2))
        self._bl_label_nome = tk.Label(self.f_baseline, text="",
                                       font=('Consolas', 9, 'italic'), bg=bg,
                                       wraplength=360, justify='left')
        self._bl_label_nome.pack(padx=10)

        ctrl = tk.Frame(self.f_baseline, bg=bg)
        ctrl.pack(fill=tk.X, padx=10, pady=6)

        def anchor_row(label, var, default):
            row = tk.Frame(ctrl, bg=bg)
            row.pack(fill=tk.X, pady=3)
            tk.Label(row, text=label, width=9, anchor='w', bg=bg, font=('Consolas', 9)).pack(side=tk.LEFT)
            var.set(default)
            e = tk.Entry(row, textvariable=var, width=8, font=('Consolas', 9))
            e.pack(side=tk.LEFT, padx=4)
            e.bind('<Return>',   lambda ev: self._bl_from_entry())
            e.bind('<FocusOut>', lambda ev: self._bl_from_entry())
            tk.Label(row, text='cm⁻¹', bg=bg, font=('Consolas', 8), fg='#555555').pack(side=tk.LEFT)

        self._bl_x1_var = tk.StringVar()
        self._bl_x2_var = tk.StringVar()
        anchor_row('Anchor 1', self._bl_x1_var, '4000')
        anchor_row('Anchor 2', self._bl_x2_var, '2000')

        self._bl_info = tk.Label(self.f_baseline, text="", font=('Consolas', 8),
                                 bg=bg, fg='#555555', justify='left')
        self._bl_info.pack(padx=10, pady=(0, 4))

        tk.Label(self.f_baseline, text="(trascina le linee o digita i valori)",
                 font=('Consolas', 8, 'italic'), bg=bg, fg='#888888').pack()

        btn_bl = tk.Frame(self.f_baseline, bg=bg)
        btn_bl.pack(fill=tk.X, padx=10, pady=10)
        tk.Button(btn_bl, text="Apply", command=self._bl_applica,
                  bg='#2ecc71', fg='white', font=('Arial', 9, 'bold')).pack(side=tk.LEFT, expand=True, fill=tk.X, padx=2)
        tk.Button(btn_bl, text="Cancel", command=self._bl_annulla).pack(side=tk.LEFT, expand=True, fill=tk.X, padx=2)

    # --- BASELINE: apertura pannello ---
    def apri_baseline(self):
        sel = self._selezione_effettiva()
        if len(sel) != 1:
            messagebox.showwarning("Baseline", "Seleziona esattamente uno spettro FTIR.")
            return
        name = self.file_listbox.get(sel[0])
        if self.spectra[name]['info']['Type'] != 'FTIR':
            messagebox.showwarning("Baseline", "La baseline lineare è pensata per spettri FTIR.")
            return

        self._bl_name = name
        self._bl_x = self.spectra[name]['df'].index.to_numpy(dtype=float)
        self._bl_y = self.spectra[name]['df'][name].to_numpy(dtype=float)
        self._bl_label_nome.config(text=name)

        # Ancore di default clampate al range reale dello spettro
        x_min, x_max = self._bl_x.min(), self._bl_x.max()
        a1 = min(4000.0, x_max)
        a2 = 2000.0 if x_min <= 2000.0 <= x_max else (x_min + x_max) / 2
        self._bl_syncing = True
        self._bl_x1_var.set(f"{a1:.0f}")
        self._bl_x2_var.set(f"{a2:.0f}")
        self._bl_syncing = False

        # Prepara grafico
        self.ax.clear()
        self.v_line = self.v_text = None
        self.ax.set_xlabel("Wavenumber (cm⁻¹)")
        self.ax.set_ylabel(self._ftir_ylabel([name]))
        self.ax.grid(True, linestyle=':', alpha=0.6)

        self._bl_lines['orig'], = self.ax.plot(self._bl_x, self._bl_y,
                                               color='steelblue', alpha=0.5, lw=1.5, label=f'{name} (orig)')
        self._bl_lines['base'], = self.ax.plot(self._bl_x, np.zeros_like(self._bl_y),
                                               color='orange', lw=1.2, ls='--', label='Baseline')
        self._bl_lines['res'],  = self.ax.plot(self._bl_x, self._bl_y.copy(),
                                               color='green', lw=2, label=f'{name}_blc')
        self.ax.legend(fontsize=8)

        # Linee d'ancoraggio trascinabili
        self._bl_va = self.ax.axvline(x=a1, color='purple', ls='--', lw=1.5, alpha=0.8)
        self._bl_vb = self.ax.axvline(x=a2, color='purple', ls='--', lw=1.5, alpha=0.8)
        self._bl_drag = None
        self._bl_cids = [
            self.canvas.mpl_connect('button_press_event',   self._bl_drag_start),
            self.canvas.mpl_connect('motion_notify_event',  self._bl_drag_move),
            self.canvas.mpl_connect('button_release_event', self._bl_drag_stop),
        ]
        self.ax.invert_xaxis()
        self._bl_aggiorna_preview(rescale=True)
        self.canvas.draw()

        self.f_metadata.pack_forget()
        self.f_baseline.pack(fill=tk.BOTH, expand=True)

    # --- BASELINE: calcolo retta passante per le due ancore ---
    def _bl_baseline_array(self):
        x, y = self._bl_x, self._bl_y
        xa = self._bl_va.get_xdata()[0]
        xb = self._bl_vb.get_xdata()[0]
        ia = int(np.argmin(np.abs(x - xa)))
        ib = int(np.argmin(np.abs(x - xb)))
        xa, ya = x[ia], y[ia]
        xb, yb = x[ib], y[ib]
        if xa == xb:
            base = np.full_like(y, ya)
        else:
            # Retta tra le due ancore; costante oltre ciascun estremo (niente
            # estrapolazione della pendenza, che distorcerebbe lo spettro)
            if xa < xb:
                base = np.interp(x, [xa, xb], [ya, yb])
            else:
                base = np.interp(x, [xb, xa], [yb, ya])
        return base, (xa, ya), (xb, yb)

    # --- BASELINE: preview live ---
    def _bl_aggiorna_preview(self, rescale=False):
        if self._bl_x is None or 'base' not in self._bl_lines:
            return
        base, (xa, ya), (xb, yb) = self._bl_baseline_array()
        corr = self._bl_y - base
        self._bl_lines['base'].set_ydata(base)
        self._bl_lines['res'].set_ydata(corr)
        self._bl_info.config(text=f"A1: {xa:.0f} = {ya:.4f}\nA2: {xb:.0f} = {yb:.4f}")
        if rescale:
            y_all  = np.concatenate([self._bl_y, corr, base])
            margin = (y_all.max() - y_all.min()) * 0.05 or 0.05
            self.ax.set_ylim(y_all.min() - margin, y_all.max() + margin)
        self.canvas.draw_idle()

    # --- BASELINE: sincronizzazione caselle → linee ---
    def _bl_from_entry(self):
        if self._bl_syncing or self._bl_x is None:
            return
        try:
            a1 = float(self._bl_x1_var.get().replace(',', '.'))
            a2 = float(self._bl_x2_var.get().replace(',', '.'))
        except ValueError:
            return
        self._bl_va.set_xdata([a1, a1])
        self._bl_vb.set_xdata([a2, a2])
        self._bl_aggiorna_preview()

    # --- BASELINE: drag delle linee d'ancoraggio ---
    def _bl_drag_start(self, event):
        if event.inaxes != self.ax or event.button != 1:
            return
        for key, line in [('a', self._bl_va), ('b', self._bl_vb)]:
            x_disp = self.ax.transData.transform((line.get_xdata()[0], 0))[0]
            if abs(event.x - x_disp) < 8:
                self._bl_drag = key
                return

    def _bl_drag_move(self, event):
        if self._bl_drag is None or event.inaxes != self.ax or event.xdata is None:
            return
        xv = event.xdata
        line = self._bl_va if self._bl_drag == 'a' else self._bl_vb
        var  = self._bl_x1_var if self._bl_drag == 'a' else self._bl_x2_var
        line.set_xdata([xv, xv])
        self._bl_syncing = True
        var.set(f"{xv:.0f}")
        self._bl_syncing = False
        self._bl_aggiorna_preview()

    def _bl_drag_stop(self, event):
        self._bl_drag = None

    # --- BASELINE: applica ---
    def _bl_applica(self):
        base, (xa, ya), (xb, yb) = self._bl_baseline_array()
        corr = self._round_sig(self._bl_y - base, sig=6)
        new_name = f"{self._bl_name}_blc"
        final_name = new_name
        counter = 2
        while final_name in self.spectra:
            final_name = f"{new_name}_{counter}"
            counter += 1
        new_info = dict(self.spectra[self._bl_name]['info'])
        new_info['Sample']          = final_name
        new_info['Derived from']    = self._bl_name
        new_info['Operation']       = 'Baseline (linear)'
        new_info['Baseline Anchor 1'] = f"{xa:.1f}"
        new_info['Baseline Anchor 2'] = f"{xb:.1f}"
        self.spectra[final_name] = {
            'df':   pd.DataFrame({final_name: corr}, index=self.spectra[self._bl_name]['df'].index),
            'info': new_info,
        }
        self._dirty = True
        self._bl_chiudi_pannello()
        self.aggiorna_vista()

    # --- BASELINE: annulla ---
    def _bl_annulla(self):
        self._bl_chiudi_pannello()
        self.aggiorna_vista()

    # --- BASELINE: chiude il pannello ---
    def _bl_chiudi_pannello(self):
        for cid in self._bl_cids:
            self.canvas.mpl_disconnect(cid)
        self._bl_cids  = []
        self._bl_drag  = None
        self._bl_va    = None
        self._bl_vb    = None
        self._bl_name  = None
        self._bl_x     = None
        self._bl_y     = None
        self._bl_lines = {}
        self.f_baseline.pack_forget()
        self.f_metadata.pack(fill=tk.BOTH, expand=True)

    # --- PANNELLO BASELINE ADATTIVA: costruzione widget ---
    def _costruisci_pannello_adattivo(self):
        bg = '#f0f4f7'

        tk.Label(self.f_adaptive, text="ADAPTIVE BASELINE",
                 font=('Arial', 10, 'bold'), bg=bg).pack(pady=(15, 2))
        self._ab_label_nome = tk.Label(self.f_adaptive, text="",
                                       font=('Consolas', 9, 'italic'), bg=bg,
                                       wraplength=360, justify='left')
        self._ab_label_nome.pack(padx=10)

        ctrl = tk.Frame(self.f_adaptive, bg=bg)
        ctrl.pack(fill=tk.X, padx=10, pady=8)
        row = tk.Frame(ctrl, bg=bg)
        row.pack(fill=tk.X, pady=4)
        tk.Label(row, text='Adapt', width=6, anchor='w', bg=bg, font=('Consolas', 9)).pack(side=tk.LEFT)
        self._ab_s = tk.DoubleVar(value=0.0)
        tk.Scale(row, variable=self._ab_s, from_=0.0, to=1.0, resolution=0.01,
                 orient=tk.HORIZONTAL, length=250, bg=bg,
                 font=('Consolas', 8), showvalue=True).pack(side=tk.LEFT)
        self._ab_s.trace_add('write', self._ab_aggiorna_preview)

        tk.Label(self.f_adaptive, text="◀ piatta          aderente ▶",
                 font=('Consolas', 8), bg=bg, fg='#555555').pack()
        self._ab_info = tk.Label(self.f_adaptive, text="", font=('Consolas', 8),
                                 bg=bg, fg='#555555', justify='left')
        self._ab_info.pack(padx=10, pady=(2, 4))

        btn = tk.Frame(self.f_adaptive, bg=bg)
        btn.pack(fill=tk.X, padx=10, pady=10)
        tk.Button(btn, text="Apply", command=self._ab_applica,
                  bg='#2ecc71', fg='white', font=('Arial', 9, 'bold')).pack(side=tk.LEFT, expand=True, fill=tk.X, padx=2)
        tk.Button(btn, text="Cancel", command=self._ab_annulla).pack(side=tk.LEFT, expand=True, fill=tk.X, padx=2)

    @staticmethod
    def _smooth_boxcar(y, window, passes=3):
        """Media mobile ripetuta (≈ gaussiana per il TLC) via cumsum, O(N).
        Padding riflesso per ridurre le distorsioni ai bordi."""
        n = len(y)
        window = int(max(1, min(window, n)))
        if window <= 1:
            return y.astype(float).copy()
        if window % 2 == 0:
            window += 1
        half = window // 2
        out = y.astype(float)
        for _ in range(passes):
            padded = np.pad(out, half, mode='reflect')
            csum = np.cumsum(np.insert(padded, 0, 0.0))
            out = (csum[window:] - csum[:-window]) / window
        return out

    @staticmethod
    def _savgol(y, window, polyorder=2):
        """Filtro Savitzky-Golay puro numpy: fit polinomiale locale, preserva
        altezza e forma dei picchi meglio della media mobile.
        I coefficienti di convoluzione sono la prima riga di (AᵀA)⁻¹Aᵀ."""
        n = len(y)
        window = int(window)
        if window % 2 == 0:
            window += 1
        if window > n:
            window = n if n % 2 == 1 else n - 1
        if window < 3 or window <= polyorder + 1:
            return y.astype(float).copy()
        half = window // 2
        j = np.arange(-half, half + 1)
        A = np.vander(j, polyorder + 1, increasing=True)   # colonne j^0 … j^p
        kernel = np.linalg.pinv(A)[0]                       # coeff. per il valore centrale
        # Riflessione "dispari": estende il trend lineare ai bordi senza gomiti
        padded = np.pad(y.astype(float), half, mode='reflect', reflect_type='odd')
        return np.convolve(padded, kernel[::-1], mode='valid')

    # --- BASELINE ADATTIVA: finestra dallo slider (log) ---
    def _ab_window(self):
        n = len(self._ab_y)
        s = self._ab_s.get()               # 0 → finestra n (piatta); 1 → finestra 1 (segue)
        return int(max(1, round(n ** (1.0 - s))))

    def _ab_baseline_array(self):
        n = len(self._ab_y)
        w = self._ab_window()
        y = self._ab_y
        if w >= n:
            # Estremo "piatta": retta orizzontale sotto lo spettro
            return np.full_like(y, y.min())
        # Vincolo riferito a una versione leggermente denoised dello spettro: così
        # la baseline sta sotto il segnale ma non insegue i minimi del rumore
        # (nelle zone piatte passa per il centro del rumore → il risultato può
        #  scendere sotto zero, che è rumore sperimentale legittimo).
        w_ref = max(1, w // 8)
        y_ref = self._smooth_boxcar(y, w_ref) if w_ref > 1 else y
        z = y_ref.copy()
        for _ in range(12):
            z = self._smooth_boxcar(z, w)
            z = np.minimum(z, y_ref)
        return z

    # --- BASELINE ADATTIVA: apertura pannello ---
    def apri_baseline_adattiva(self):
        sel = self._selezione_effettiva()
        if len(sel) != 1:
            messagebox.showwarning("Adaptive Baseline", "Seleziona esattamente uno spettro.")
            return
        name = self.file_listbox.get(sel[0])
        df = self.spectra[name]['df']
        self._ab_name = name
        self._ab_x = df.index.to_numpy(dtype=float)
        self._ab_y = df[name].to_numpy(dtype=float)
        self._ab_label_nome.config(text=name)
        self._ab_s.set(0.0)

        self.ax.clear()
        self.v_line = self.v_text = None
        tipo = self.spectra[name]['info']['Type']
        if tipo == 'FTIR':
            self.ax.set_xlabel("Wavenumber (cm⁻¹)")
            self.ax.set_ylabel(self._ftir_ylabel([name]))
        elif tipo == 'Raman':
            self.ax.set_xlabel("Raman shift (cm⁻¹)")
            self.ax.set_ylabel("Intensity (A.U.)")
        elif tipo == 'Fluorescence':
            self.ax.set_xlabel("Wavelength (nm)")
            self.ax.set_ylabel("Intensity (A.U.)")
        else:
            self.ax.set_xlabel("Wavelength (nm)")
            self.ax.set_ylabel("Absorbance (A)")
        self.ax.grid(True, linestyle=':', alpha=0.6)

        # In anteprima mostriamo solo originale + baseline (stessa scala):
        # il risultato corretto appare solo dopo Apply, per non schiacciare l'asse Y
        self._ab_lines['orig'], = self.ax.plot(self._ab_x, self._ab_y,
                                               color='steelblue', alpha=0.6, lw=1.5, label=f'{name} (orig)')
        self._ab_lines['base'], = self.ax.plot(self._ab_x, np.zeros_like(self._ab_y),
                                               color='orange', lw=1.5, ls='--', label='Baseline')
        self.ax.legend(fontsize=8)
        if tipo == 'FTIR':
            self.ax.invert_xaxis()
        self._ab_aggiorna_preview()
        # Fissa la scala Y su quella dello spettro originale (la baseline vi rientra)
        margin = (self._ab_y.max() - self._ab_y.min()) * 0.05 or 0.05
        self.ax.set_ylim(self._ab_y.min() - margin, self._ab_y.max() + margin)
        self.canvas.draw()

        self.f_metadata.pack_forget()
        self.f_adaptive.pack(fill=tk.BOTH, expand=True)

    # --- BASELINE ADATTIVA: preview live (solo originale + baseline) ---
    def _ab_aggiorna_preview(self, *_):
        if self._ab_x is None or 'base' not in self._ab_lines:
            return
        base = self._ab_baseline_array()
        self._ab_lines['base'].set_ydata(base)
        self._ab_info.config(text=f"finestra ≈ {self._ab_window()} punti")
        self.canvas.draw_idle()

    # --- BASELINE ADATTIVA: applica ---
    def _ab_applica(self):
        base = self._ab_baseline_array()
        corr = self._round_sig(self._ab_y - base, sig=6)
        new_name = f"{self._ab_name}_abl"
        final_name = new_name
        counter = 2
        while final_name in self.spectra:
            final_name = f"{new_name}_{counter}"
            counter += 1
        new_info = dict(self.spectra[self._ab_name]['info'])
        new_info['Sample']         = final_name
        new_info['Derived from']   = self._ab_name
        new_info['Operation']      = 'Adaptive Baseline'
        new_info['Adapt']          = f"{self._ab_s.get():.2f}"
        new_info['Baseline Window'] = f"{self._ab_window()} pts"
        self.spectra[final_name] = {
            'df':   pd.DataFrame({final_name: corr}, index=self.spectra[self._ab_name]['df'].index),
            'info': new_info,
        }
        self._dirty = True
        self._ab_chiudi_pannello()
        self.aggiorna_vista()

    # --- BASELINE ADATTIVA: annulla ---
    def _ab_annulla(self):
        self._ab_chiudi_pannello()
        self.aggiorna_vista()

    # --- BASELINE ADATTIVA: chiude il pannello ---
    def _ab_chiudi_pannello(self):
        self._ab_name  = None
        self._ab_x     = None
        self._ab_y     = None
        self._ab_lines = {}
        self.f_adaptive.pack_forget()
        self.f_metadata.pack(fill=tk.BOTH, expand=True)

    # --- PANNELLO SMOOTHING: costruzione widget ---
    def _costruisci_pannello_smoothing(self):
        bg = '#f0f4f7'

        tk.Label(self.f_smooth, text="SMOOTHING",
                 font=('Arial', 10, 'bold'), bg=bg).pack(pady=(15, 2))
        self._sm_label_nome = tk.Label(self.f_smooth, text="",
                                       font=('Consolas', 9, 'italic'), bg=bg,
                                       wraplength=360, justify='left')
        self._sm_label_nome.pack(padx=10)

        ctrl = tk.Frame(self.f_smooth, bg=bg)
        ctrl.pack(fill=tk.X, padx=10, pady=8)
        row = tk.Frame(ctrl, bg=bg)
        row.pack(fill=tk.X, pady=4)
        tk.Label(row, text='Amount', width=6, anchor='w', bg=bg, font=('Consolas', 9)).pack(side=tk.LEFT)
        self._sm_s = tk.DoubleVar(value=0.0)
        tk.Scale(row, variable=self._sm_s, from_=0.0, to=1.0, resolution=0.01,
                 orient=tk.HORIZONTAL, length=250, bg=bg,
                 font=('Consolas', 8), showvalue=True).pack(side=tk.LEFT)
        self._sm_s.trace_add('write', self._sm_aggiorna_preview)

        tk.Label(self.f_smooth, text="◀ nessuno          forte ▶",
                 font=('Consolas', 8), bg=bg, fg='#555555').pack()

        # Metodo di smoothing
        met_row = tk.Frame(ctrl, bg=bg)
        met_row.pack(fill=tk.X, pady=(8, 2))
        tk.Label(met_row, text='Method', width=6, anchor='w', bg=bg, font=('Consolas', 9)).pack(side=tk.LEFT)
        self._sm_method = tk.StringVar(value='movavg')
        tk.Radiobutton(met_row, text='Moving avg', variable=self._sm_method, value='movavg',
                       bg=bg, font=('Consolas', 8), command=self._sm_aggiorna_preview).pack(side=tk.LEFT)
        tk.Radiobutton(met_row, text='Sav-Golay', variable=self._sm_method, value='savgol',
                       bg=bg, font=('Consolas', 8), command=self._sm_aggiorna_preview).pack(side=tk.LEFT)

        # Ordine del polinomio (solo Savitzky-Golay)
        poly_row = tk.Frame(ctrl, bg=bg)
        poly_row.pack(fill=tk.X, pady=2)
        tk.Label(poly_row, text='Poly', width=6, anchor='w', bg=bg, font=('Consolas', 9)).pack(side=tk.LEFT)
        self._sm_poly = tk.IntVar(value=2)
        for p_val in (2, 3, 4):
            tk.Radiobutton(poly_row, text=str(p_val), variable=self._sm_poly, value=p_val,
                           bg=bg, font=('Consolas', 8), command=self._sm_aggiorna_preview).pack(side=tk.LEFT, padx=4)

        self._sm_info = tk.Label(self.f_smooth, text="", font=('Consolas', 8),
                                 bg=bg, fg='#555555', justify='left')
        self._sm_info.pack(padx=10, pady=(2, 4))

        btn = tk.Frame(self.f_smooth, bg=bg)
        btn.pack(fill=tk.X, padx=10, pady=10)
        tk.Button(btn, text="Apply", command=self._sm_applica,
                  bg='#2ecc71', fg='white', font=('Arial', 9, 'bold')).pack(side=tk.LEFT, expand=True, fill=tk.X, padx=2)
        tk.Button(btn, text="Cancel", command=self._sm_annulla).pack(side=tk.LEFT, expand=True, fill=tk.X, padx=2)

    # --- SMOOTHING: finestra dallo slider ---
    def _sm_window(self):
        n = len(self._sm_y)
        s = self._sm_s.get()
        maxw = max(5, n // 10)          # finestra massima ~10% dei punti
        return int(max(1, round(maxw ** s)))   # s=0 → 1 (nessuno); s=1 → maxw

    def _sm_smoothed(self):
        w = self._sm_window()
        if self._sm_method.get() == 'savgol':
            p = self._sm_poly.get()
            win = max(w, p + 2)
            if win % 2 == 0:
                win += 1
            return self._savgol(self._sm_y, win, p)
        if w <= 1:
            return self._sm_y.copy()
        return self._smooth_boxcar(self._sm_y, w)

    def _sm_window_effettiva(self):
        """Finestra realmente usata (per l'etichetta e i metadati)."""
        w = self._sm_window()
        if self._sm_method.get() == 'savgol':
            p = self._sm_poly.get()
            win = max(w, p + 2)
            if win % 2 == 0:
                win += 1
            return win
        return w

    # --- SMOOTHING: apertura pannello ---
    def apri_smoothing(self):
        sel = self._selezione_effettiva()
        if len(sel) != 1:
            messagebox.showwarning("Smoothing", "Seleziona esattamente uno spettro.")
            return
        name = self.file_listbox.get(sel[0])
        df = self.spectra[name]['df']
        self._sm_name = name
        self._sm_x = df.index.to_numpy(dtype=float)
        self._sm_y = df[name].to_numpy(dtype=float)
        self._sm_label_nome.config(text=name)
        self._sm_s.set(0.0)

        self.ax.clear()
        self.v_line = self.v_text = None
        tipo = self.spectra[name]['info']['Type']
        if tipo == 'FTIR':
            self.ax.set_xlabel("Wavenumber (cm⁻¹)")
            self.ax.set_ylabel(self._ftir_ylabel([name]))
        elif tipo == 'Raman':
            self.ax.set_xlabel("Raman shift (cm⁻¹)")
            self.ax.set_ylabel("Intensity (A.U.)")
        elif tipo == 'Fluorescence':
            self.ax.set_xlabel("Wavelength (nm)")
            self.ax.set_ylabel("Intensity (A.U.)")
        else:
            self.ax.set_xlabel("Wavelength (nm)")
            self.ax.set_ylabel("Absorbance (A)")
        self.ax.grid(True, linestyle=':', alpha=0.6)

        self._sm_lines['orig'], = self.ax.plot(self._sm_x, self._sm_y,
                                               color='steelblue', alpha=0.4, lw=1.2, label=f'{name} (orig)')
        self._sm_lines['smooth'], = self.ax.plot(self._sm_x, self._sm_y.copy(),
                                                 color='green', lw=1.8, label=f'{name}_smooth')
        self.ax.legend(fontsize=8)
        if tipo == 'FTIR':
            self.ax.invert_xaxis()
        self._sm_aggiorna_preview()
        # Scala Y fissa su quella dell'originale (lo smoothing vi rientra)
        margin = (self._sm_y.max() - self._sm_y.min()) * 0.05 or 0.05
        self.ax.set_ylim(self._sm_y.min() - margin, self._sm_y.max() + margin)
        self.canvas.draw()

        self.f_metadata.pack_forget()
        self.f_smooth.pack(fill=tk.BOTH, expand=True)

    # --- SMOOTHING: preview live ---
    def _sm_aggiorna_preview(self, *_):
        if self._sm_x is None or 'smooth' not in self._sm_lines:
            return
        self._sm_lines['smooth'].set_ydata(self._sm_smoothed())
        if self._sm_method.get() == 'savgol':
            self._sm_info.config(text=f"Sav-Golay: finestra {self._sm_window_effettiva()} pt, ordine {self._sm_poly.get()}")
        else:
            self._sm_info.config(text=f"Moving avg: finestra ≈ {self._sm_window_effettiva()} pt")
        self.canvas.draw_idle()

    # --- SMOOTHING: applica ---
    def _sm_applica(self):
        smoothed = self._round_sig(self._sm_smoothed(), sig=6)
        new_name = f"{self._sm_name}_smooth"
        final_name = new_name
        counter = 2
        while final_name in self.spectra:
            final_name = f"{new_name}_{counter}"
            counter += 1
        metodo = ('Savitzky-Golay' if self._sm_method.get() == 'savgol' else 'Moving average')
        new_info = dict(self.spectra[self._sm_name]['info'])
        new_info['Sample']         = final_name
        new_info['Derived from']   = self._sm_name
        new_info['Operation']      = 'Smoothing'
        new_info['Smooth Method']  = metodo
        new_info['Smooth Window']  = f"{self._sm_window_effettiva()} pts"
        if self._sm_method.get() == 'savgol':
            new_info['Smooth Poly'] = str(self._sm_poly.get())
        self.spectra[final_name] = {
            'df':   pd.DataFrame({final_name: smoothed}, index=self.spectra[self._sm_name]['df'].index),
            'info': new_info,
        }
        self._dirty = True
        self._sm_chiudi_pannello()
        self.aggiorna_vista()

    # --- SMOOTHING: annulla ---
    def _sm_annulla(self):
        self._sm_chiudi_pannello()
        self.aggiorna_vista()

    # --- SMOOTHING: chiude il pannello ---
    def _sm_chiudi_pannello(self):
        self._sm_name  = None
        self._sm_x     = None
        self._sm_y     = None
        self._sm_lines = {}
        self.f_smooth.pack_forget()
        self.f_metadata.pack(fill=tk.BOTH, expand=True)

    # --- PANNELLO DECONVOLUZIONE: costruzione widget ---
    @staticmethod
    def _lorentzian(x, A, x0, w):
        """Lorentziana: A = altezza del picco, x0 = centro, w = FWHM."""
        hw = w / 2.0
        return A * hw * hw / ((x - x0) ** 2 + hw * hw)

    @staticmethod
    def _pseudovoigt(x, A, x0, w, eta):
        """pseudo-Voigt (altezza-normalizzata): A[η·Lorentz + (1−η)·Gauss].
        η=1 → lorentziana pura; η=0 → gaussiana pura. Stessa FWHM w, stessa altezza A."""
        hw = w / 2.0
        lor = hw * hw / ((x - x0) ** 2 + hw * hw)
        gau = np.exp(-4.0 * np.log(2.0) * (x - x0) ** 2 / (w * w))
        return A * (eta * lor + (1.0 - eta) * gau)

    @staticmethod
    def _pv_area(A, w, eta):
        """Area analitica di una pseudo-Voigt altezza-normalizzata."""
        lor = np.pi * w / 2.0
        gau = 0.5 * w * np.sqrt(np.pi / np.log(2.0))
        return A * (eta * lor + (1.0 - eta) * gau)

    def _costruisci_pannello_deconv(self):
        bg = '#f0f4f7'

        tk.Label(self.f_deconv, text="DECONVOLUTION (Lorentzian)",
                 font=('Arial', 10, 'bold'), bg=bg).pack(pady=(15, 2))
        self._dc_label_nome = tk.Label(self.f_deconv, text="",
                                       font=('Consolas', 9, 'italic'), bg=bg,
                                       wraplength=360, justify='left')
        self._dc_label_nome.pack(padx=10)

        # Modalità aggiunta/rimozione picchi: un tasto dedicato invece di un click
        # sempre "armato", così si può guardare/zoomare il grafico senza piazzare
        # picchi per sbaglio, ed è un punto naturale per sbloccare pan/zoom della
        # toolbar se erano rimasti attivi (altrimenti il click risulta silenziosamente
        # ignorato e sembra che "non succeda nulla").
        self._dc_add_mode = tk.BooleanVar(value=True)
        tk.Checkbutton(self.f_deconv, text="✏ Add/remove peaks (click on the graph)",
                       variable=self._dc_add_mode, bg=bg, font=('Consolas', 9),
                       command=self._dc_disattiva_toolbar_mode).pack(padx=10, pady=(4, 0), anchor='w')
        tk.Label(self.f_deconv,
                 text="Left click: add peak   •   Right click: remove nearest",
                 font=('Consolas', 8), bg=bg, fg='#555555', justify='left').pack(padx=10, pady=(0, 2))

        # Zona di interesse per il fit: due linee trascinabili sul grafico,
        # oppure inserimento diretto dei valori (stesso pattern del pannello Trim)
        tk.Label(self.f_deconv, text="Fit range:",
                 font=('Consolas', 8), bg=bg, fg='#555555', justify='left').pack(padx=10, anchor='w')
        dc_range_row = tk.Frame(self.f_deconv, bg=bg)
        dc_range_row.pack(padx=10, pady=(2, 0))
        tk.Label(dc_range_row, text="From", bg=bg, font=('Consolas', 9)).pack(side=tk.LEFT)
        self._dc_entry_lo = tk.Entry(dc_range_row, width=9, font=('Consolas', 9), justify='center')
        self._dc_entry_lo.pack(side=tk.LEFT, padx=(4, 12))
        tk.Label(dc_range_row, text="To", bg=bg, font=('Consolas', 9)).pack(side=tk.LEFT)
        self._dc_entry_hi = tk.Entry(dc_range_row, width=9, font=('Consolas', 9), justify='center')
        self._dc_entry_hi.pack(side=tk.LEFT, padx=(4, 0))
        for entry in (self._dc_entry_lo, self._dc_entry_hi):
            entry.bind('<Return>', self._dc_entry_commit)
            entry.bind('<FocusOut>', self._dc_entry_commit)
        self._dc_range_label = tk.Label(self.f_deconv, text="—",
                                        font=('Consolas', 8), bg=bg, fg='#555555')
        self._dc_range_label.pack(pady=(2, 0))

        # Profilo di forma
        shape_row = tk.Frame(self.f_deconv, bg=bg)
        shape_row.pack(fill=tk.X, padx=10, pady=(2, 0))
        tk.Label(shape_row, text='Profile', width=7, anchor='w', bg=bg, font=('Consolas', 9)).pack(side=tk.LEFT)
        self._dc_shape = tk.StringVar(value='pv')
        for lbl, val in [('Lorentz', 'lorentz'), ('Gauss', 'gauss'), ('p-Voigt', 'pv')]:
            tk.Radiobutton(shape_row, text=lbl, variable=self._dc_shape, value=val,
                           bg=bg, font=('Consolas', 8),
                           command=self._dc_cambia_shape).pack(side=tk.LEFT)

        btn1 = tk.Frame(self.f_deconv, bg=bg)
        btn1.pack(fill=tk.X, padx=10, pady=(2, 0))
        tk.Button(btn1, text="Fit", command=self._dc_fit,
                  bg='#42a5f5', fg='white', font=('Arial', 9, 'bold')).pack(side=tk.LEFT, expand=True, fill=tk.X, padx=2)
        tk.Button(btn1, text="Clear peaks", command=self._dc_clear).pack(side=tk.LEFT, expand=True, fill=tk.X, padx=2)

        # Risultati: Text selezionabile + copia in clipboard
        self._dc_info = tk.Text(self.f_deconv, height=8, font=('Consolas', 8),
                                bg='#ffffff', fg='#333333', wrap='none', bd=1, relief='solid')
        self._dc_info.pack(padx=10, pady=(6, 2), fill=tk.X)
        self._dc_set_info("Nessun picco.")
        tk.Button(self.f_deconv, text="Copy results", command=self._dc_copia_risultati,
                  font=('Arial', 8)).pack(padx=10, pady=(0, 4), anchor='e')

        btn2 = tk.Frame(self.f_deconv, bg=bg)
        btn2.pack(fill=tk.X, padx=10, pady=10)
        tk.Button(btn2, text="Apply", command=self._dc_applica,
                  bg='#2ecc71', fg='white', font=('Arial', 9, 'bold')).pack(side=tk.LEFT, expand=True, fill=tk.X, padx=2)
        tk.Button(btn2, text="Cancel", command=self._dc_annulla).pack(side=tk.LEFT, expand=True, fill=tk.X, padx=2)

    # --- DECONVOLUZIONE: apertura pannello ---
    def apri_deconvoluzione(self):
        sel = self._selezione_effettiva()
        if len(sel) != 1:
            messagebox.showwarning("Deconvolution", "Seleziona esattamente uno spettro.")
            return
        name = self.file_listbox.get(sel[0])
        df = self.spectra[name]['df']
        self._dc_name = name
        self._dc_x = df.index.to_numpy(dtype=float)
        self._dc_y = df[name].to_numpy(dtype=float)
        self._dc_peaks  = []
        self._dc_offset = 0.0
        self._dc_fitted = False
        self._dc_range  = [float(self._dc_x.min()), float(self._dc_x.max())]
        self._dc_drag   = None
        self._dc_drag_moved = False
        self._dc_add_mode.set(True)
        self._dc_label_nome.config(text=name)
        self._dc_disattiva_toolbar_mode()   # pan/zoom rimasto attivo da un altro pannello bloccherebbe i click

        self._dc_cids = [
            self.canvas.mpl_connect('button_press_event',   self._dc_press),
            self.canvas.mpl_connect('motion_notify_event',  self._dc_drag_move),
            self.canvas.mpl_connect('button_release_event', self._dc_drag_stop),
        ]
        self._dc_ridisegna()

        self.f_metadata.pack_forget()
        self.f_deconv.pack(fill=tk.BOTH, expand=True)

    # --- DECONVOLUZIONE: gestione mouse (linee intervallo + picchi) ---
    def _dc_disattiva_toolbar_mode(self):
        """Spegne pan/zoom della toolbar se attivi: altrimenti il primo click sul grafico
        viene ignorato in silenzio (matplotlib lo intercetta per il pan/zoom) e sembra che
        il click 'non faccia nulla' — capita spesso dopo aver usato lo zoom per ispezionare
        lo spettro prima di piazzare i picchi."""
        mode = getattr(self.toolbar.mode, 'value', self.toolbar.mode)
        if mode == 'pan/zoom':
            self.toolbar.pan()
        elif mode == 'zoom rect':
            self.toolbar.zoom()

    def _dc_aggiungi_picco(self, x0):
        x0 = float(x0)
        idx = int(np.argmin(np.abs(self._dc_x - x0)))
        A = max(float(self._dc_y[idx]) - self._dc_offset, 1e-6)
        w = (self._dc_x.max() - self._dc_x.min()) * 0.02
        self._dc_peaks.append({'x0': x0, 'A': A, 'w': w, 'eta': self._dc_eta_default()})
        self._dc_fitted = False
        self._dc_ridisegna()

    def _dc_press(self, event):
        if event.inaxes != self.ax or event.xdata is None:
            return
        if self.toolbar.mode:        # ignora se zoom/pan attivo
            return
        # Priorità: se vicino a una linea dell'intervallo, inizia il trascinamento (ma se il
        # mouse non si sposta prima del rilascio, _dc_drag_stop lo tratta come un click e
        # aggiunge comunque un picco lì — altrimenti un peak vicino al bordo del range,
        # ora esteso di default a tutto lo spettro, non sarebbe mai selezionabile)
        if event.button == 1:
            for key, line in [('a', self._dc_va), ('b', self._dc_vb)]:
                if line is not None:
                    x_disp = self.ax.transData.transform((line.get_xdata()[0], 0))[0]
                    if abs(event.x - x_disp) < 8:
                        self._dc_drag = key
                        self._dc_drag_moved = False
                        return
            # altrimenti aggiungi un picco (se la modalità è attiva)
            if self._dc_add_mode.get():
                self._dc_aggiungi_picco(event.xdata)
        elif event.button == 3 and self._dc_peaks and self._dc_add_mode.get():
            centers = np.array([pk['x0'] for pk in self._dc_peaks])
            i = int(np.argmin(np.abs(centers - event.xdata)))
            del self._dc_peaks[i]
            self._dc_fitted = False
            self._dc_ridisegna()

    def _dc_drag_move(self, event):
        if self._dc_drag is None or event.inaxes != self.ax or event.xdata is None:
            return
        xv = self._dc_snap(event.xdata)
        lo, hi = min(self._dc_range), max(self._dc_range)
        if self._dc_drag == 'a':
            if xv >= hi:
                return
            self._dc_range[0] = xv
            self._dc_va.set_xdata([xv, xv])
        else:
            if xv <= lo:
                return
            self._dc_range[1] = xv
            self._dc_vb.set_xdata([xv, xv])
        self._dc_drag_moved = True
        self._dc_aggiorna_range_label()
        self.canvas.draw_idle()

    def _dc_drag_stop(self, event):
        # Click fermo vicino a una linea (nessun trascinamento effettivo): trattalo come
        # un normale click per aggiungere un picco, invece di non fare nulla.
        if (self._dc_drag is not None and not self._dc_drag_moved and self._dc_add_mode.get()
                and event is not None and event.inaxes == self.ax and event.xdata is not None):
            self._dc_aggiungi_picco(event.xdata)
        self._dc_drag = None
        self._dc_drag_moved = False

    def _dc_snap(self, xv):
        """Aggancia un valore x al punto dati più vicino (come nel pannello Trim):
        la zona di interesse per il fit cade sempre su un campione reale."""
        idx = int(np.argmin(np.abs(self._dc_x - xv)))
        return float(self._dc_x[idx])

    def _dc_aggiorna_range_label(self):
        if self._dc_range is None:
            return
        lo, hi = min(self._dc_range), max(self._dc_range)
        for entry, val in ((self._dc_entry_lo, lo), (self._dc_entry_hi, hi)):
            entry.delete(0, tk.END)
            entry.insert(0, f"{val:.4g}")
        n = int(np.sum((self._dc_x >= lo) & (self._dc_x <= hi)))
        self._dc_range_label.config(text=f"{n} points in fit range")

    # --- DECONVOLUZIONE: inserimento diretto della zona di interesse da tastiera ---
    def _dc_entry_commit(self, event=None):
        if self._dc_range is None:
            return
        try:
            lo_val = float(self._dc_entry_lo.get())
            hi_val = float(self._dc_entry_hi.get())
        except ValueError:
            self._dc_aggiorna_range_label()   # testo non numerico: ripristina i valori correnti
            return
        lo_val, hi_val = self._dc_snap(min(lo_val, hi_val)), self._dc_snap(max(lo_val, hi_val))
        if lo_val == hi_val:
            idx = int(np.argmin(np.abs(self._dc_x - lo_val)))
            if idx < len(self._dc_x) - 1:
                hi_val = float(self._dc_x[idx + 1])
            elif idx > 0:
                lo_val = float(self._dc_x[idx - 1])
        self._dc_range = [lo_val, hi_val]
        self._dc_ridisegna()

    def _dc_eta_default(self):
        return {'lorentz': 1.0, 'gauss': 0.0, 'pv': 0.5}[self._dc_shape.get()]

    def _dc_set_info(self, txt):
        self._dc_info.delete('1.0', 'end')
        self._dc_info.insert('1.0', txt)

    def _dc_copia_risultati(self):
        txt = self._dc_info.get('1.0', 'end').strip()
        if txt:
            self.root.clipboard_clear()
            self.root.clipboard_append(txt)

    def _dc_cambia_shape(self):
        # Reimposta η dei picchi al default della forma scelta e invalida il fit
        eta = self._dc_eta_default()
        for pk in self._dc_peaks:
            pk['eta'] = eta
        self._dc_fitted = False
        self._dc_ridisegna()

    # --- DECONVOLUZIONE: fit non-lineare (scipy) ---
    def _dc_fit(self):
        if not self._dc_peaks:
            messagebox.showwarning("Deconvolution", "Aggiungi almeno un picco (clic sul grafico).")
            return
        try:
            from scipy.optimize import curve_fit
        except ImportError:
            messagebox.showerror("Deconvolution",
                "Il fit richiede scipy, che non risulta installato.\n"
                "Installa con:  pip install scipy")
            return

        # Restringi il fit all'intervallo selezionato
        lo_r, hi_r = min(self._dc_range), max(self._dc_range)
        mask = (self._dc_x >= lo_r) & (self._dc_x <= hi_r)
        shape = self._dc_shape.get()
        fit_eta = (shape == 'pv')
        per = 4 if fit_eta else 3     # parametri per picco
        eta_fisso = 1.0 if shape == 'lorentz' else 0.0
        if mask.sum() < per * len(self._dc_peaks) + 1:
            messagebox.showwarning("Deconvolution",
                "L'intervallo selezionato ha troppi pochi punti per il numero di picchi.")
            return
        x, y = self._dc_x[mask], self._dc_y[mask]
        npk = len(self._dc_peaks)
        xr = x.max() - x.min()
        p0, lo, hi = [], [], []
        for pk in self._dc_peaks:
            A0  = max(pk['A'], 0.0)
            x00 = min(max(pk['x0'], x.min()), x.max())
            w0  = min(max(pk['w'], xr * 1e-4), xr)
            p0 += [A0, x00, w0]
            lo += [0.0,    x.min(), xr * 1e-4]
            hi += [np.inf, x.max(), xr]
            if fit_eta:
                e0 = min(max(pk.get('eta', 0.5), 0.0), 1.0)
                p0.append(e0); lo.append(0.0); hi.append(1.0)
        p0.append(max(self._dc_offset, 0.0)); lo.append(0.0); hi.append(np.inf)  # offset ≥ 0

        def model(xx, *p):
            s = np.full_like(xx, p[-1], dtype=float)
            for i in range(npk):
                b = i * per
                A, x0, w = p[b:b+3]
                eta = p[b+3] if fit_eta else eta_fisso
                s += self._pseudovoigt(xx, A, x0, w, eta)
            return s

        try:
            popt, _ = curve_fit(model, x, y, p0=p0, bounds=(lo, hi), maxfev=20000)
        except Exception as e:
            messagebox.showerror("Deconvolution", f"Fit non riuscito: {e}")
            return

        nuovi = []
        for i in range(npk):
            b = i * per
            eta = popt[b+3] if fit_eta else eta_fisso
            nuovi.append({'A': popt[b], 'x0': popt[b+1], 'w': abs(popt[b+2]), 'eta': eta})
        self._dc_peaks = nuovi
        self._dc_offset = float(popt[-1])
        self._dc_fitted = True
        self._dc_ridisegna()

    # --- DECONVOLUZIONE: (ri)disegna spettro, componenti e somma ---
    def _dc_ridisegna(self):
        self.ax.clear()
        self.v_line = self.v_text = None
        x, y = self._dc_x, self._dc_y
        tipo = self.spectra[self._dc_name]['info']['Type']
        if tipo == 'FTIR':
            self.ax.set_xlabel("Wavenumber (cm⁻¹)")
            self.ax.set_ylabel(self._ftir_ylabel([self._dc_name]))
        elif tipo == 'Raman':
            self.ax.set_xlabel("Raman shift (cm⁻¹)")
            self.ax.set_ylabel("Intensity (A.U.)")
        elif tipo == 'Fluorescence':
            self.ax.set_xlabel("Wavelength (nm)")
            self.ax.set_ylabel("Intensity (A.U.)")
        else:
            self.ax.set_xlabel("Wavelength (nm)")
            self.ax.set_ylabel("Absorbance (A)")
        self.ax.grid(True, linestyle=':', alpha=0.6)

        self.ax.plot(x, y, color='steelblue', alpha=0.5, lw=1.5, label=f'{self._dc_name}')

        # Il fit (somma e componenti) viene disegnato solo dentro la zona di interesse:
        # è lì che curve_fit ottimizza i parametri, fuori da quella finestra la curva non
        # ha significato ed estenderla a tutto lo spettro sarebbe fuorviante.
        lo_r, hi_r = min(self._dc_range), max(self._dc_range)
        mask = (x >= lo_r) & (x <= hi_r)

        total = np.full_like(x, self._dc_offset, dtype=float)
        colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
        righe = []
        for i, pk in enumerate(self._dc_peaks):
            eta = pk.get('eta', 1.0)
            comp = self._pseudovoigt(x, pk['A'], pk['x0'], pk['w'], eta)
            total = total + comp
            # Componente = profilo puro (sempre ≥ 0), disegnato da zero
            self.ax.plot(x[mask], comp[mask], color=colors[i % len(colors)],
                         lw=1.0, ls='--', alpha=0.8)
            self.ax.axvline(pk['x0'], color=colors[i % len(colors)], lw=0.6, alpha=0.4)
            area = self._pv_area(pk['A'], pk['w'], eta)
            eta_txt = f"  η={eta:.2f}" if self._dc_shape.get() == 'pv' else ""
            righe.append(f"P{i+1}: x0={pk['x0']:.1f}  FWHM={pk['w']:.1f}  A={pk['A']:.3f}  area={area:.2f}{eta_txt}")

        if self._dc_peaks:
            self.ax.plot(x[mask], total[mask], color='red', lw=1.8,
                         label='fit' if self._dc_fitted else 'guess')

        if self._dc_peaks and self._dc_fitted:
            ss_res = np.sum((y[mask] - total[mask]) ** 2)
            ss_tot = np.sum((y[mask] - y[mask].mean()) ** 2)
            r2 = 1 - ss_res / ss_tot if ss_tot else 0.0
            righe.append(f"offset={self._dc_offset:.3f}   R²={r2:.4f}")

        # Linee dell'intervallo di fit (trascinabili, o inseribile da tastiera)
        lo_r, hi_r = self._dc_range
        self._dc_va = self.ax.axvline(lo_r, color='purple', ls='--', lw=1.5, alpha=0.8)
        self._dc_vb = self.ax.axvline(hi_r, color='purple', ls='--', lw=1.5, alpha=0.8)
        self._dc_aggiorna_range_label()

        self.ax.legend(fontsize=8)
        if tipo == 'FTIR':
            self.ax.invert_xaxis()
        self.canvas.draw()
        righe.append(f"range: {min(self._dc_range):.0f} – {max(self._dc_range):.0f}")

        stato = "" if self._dc_fitted else "  (guess — premi Fit)"
        self._dc_set_info(('\n'.join(righe) + stato) if righe else "Nessun picco.")

    # --- DECONVOLUZIONE: azzera i picchi ---
    def _dc_clear(self):
        self._dc_peaks = []
        self._dc_offset = 0.0
        self._dc_fitted = False
        self._dc_ridisegna()

    # --- DECONVOLUZIONE: applica (salva fit + componenti) ---
    def _dc_applica(self):
        if not self._dc_peaks:
            messagebox.showwarning("Deconvolution", "Nessun picco da salvare.")
            return
        if not self._dc_fitted:
            if not messagebox.askyesno("Deconvolution",
                    "Il fit non è stato eseguito: salvo comunque i valori attuali (guess)?"):
                return
        x = self._dc_x
        idx = self.spectra[self._dc_name]['df'].index
        tipo = self.spectra[self._dc_name]['info']['Type']

        def salva(suffix, y_arr, extra):
            base = f"{self._dc_name}_{suffix}"
            final = base
            c = 2
            while final in self.spectra:
                final = f"{base}_{c}"; c += 1
            shape_nome = {'lorentz': 'Lorentzian', 'gauss': 'Gaussian',
                          'pv': 'pseudo-Voigt'}[self._dc_shape.get()]
            info = {'Sample': final, 'Date': 'N/A', 'Type': tipo,
                    'Operation': f'Deconvolution ({shape_nome})', 'Derived from': self._dc_name}
            info.update(extra)
            self.spectra[final] = {
                'df': pd.DataFrame({final: self._round_sig(y_arr, sig=6)}, index=idx),
                'info': info,
            }

        # Somma (fit totale) + ciascuna componente
        total = np.full_like(x, self._dc_offset, dtype=float)
        for pk in self._dc_peaks:
            total = total + self._pseudovoigt(x, pk['A'], pk['x0'], pk['w'], pk.get('eta', 1.0))
        salva('fit', total, {'Components': str(len(self._dc_peaks)),
                             'Offset': f"{self._dc_offset:.4f}"})
        for i, pk in enumerate(self._dc_peaks):
            eta = pk.get('eta', 1.0)
            comp = self._pseudovoigt(x, pk['A'], pk['x0'], pk['w'], eta)   # profilo puro (≥ 0)
            extra = {
                'Center':    f"{pk['x0']:.2f}",
                'FWHM':      f"{pk['w']:.2f}",
                'Amplitude': f"{pk['A']:.4f}",
                'Area':      f"{self._pv_area(pk['A'], pk['w'], eta):.3f}",
            }
            if self._dc_shape.get() == 'pv':
                extra['eta (L/G)'] = f"{eta:.3f}"
            salva(f'L{i+1}', comp, extra)

        self._dirty = True
        self._dc_chiudi_pannello()
        self.aggiorna_vista()

    # --- DECONVOLUZIONE: annulla ---
    def _dc_annulla(self):
        self._dc_chiudi_pannello()
        self.aggiorna_vista()

    # --- DECONVOLUZIONE: chiude il pannello ---
    def _dc_chiudi_pannello(self):
        for cid in self._dc_cids:
            self.canvas.mpl_disconnect(cid)
        self._dc_cids   = []
        self._dc_drag   = None
        self._dc_va     = None
        self._dc_vb     = None
        self._dc_range  = None
        self._dc_name   = None
        self._dc_x      = None
        self._dc_y      = None
        self._dc_peaks  = []
        self._dc_fitted = False
        self.f_deconv.pack_forget()
        self.f_metadata.pack(fill=tk.BOTH, expand=True)

    # --- PANNELLO TRIM: costruzione widget (una sola volta) ---
    def _costruisci_pannello_trim(self):
        bg = '#f0f4f7'

        tk.Label(self.f_trim, text="TRIM SPECTRUM",
                 font=('Arial', 10, 'bold'), bg=bg).pack(pady=(15, 2))
        self._tr_label_nome = tk.Label(self.f_trim, text="", font=('Consolas', 9, 'italic'), bg=bg)
        self._tr_label_nome.pack()

        tk.Label(self.f_trim,
                 text="Drag the two vertical lines, or type the\nstart/end values directly (Enter to confirm).",
                 font=('Consolas', 8), bg=bg, fg='#555555', justify='left').pack(padx=10, pady=(6, 2))

        range_row = tk.Frame(self.f_trim, bg=bg)
        range_row.pack(padx=10, pady=(4, 0))
        tk.Label(range_row, text="From", bg=bg, font=('Consolas', 9)).pack(side=tk.LEFT)
        self._tr_entry_lo = tk.Entry(range_row, width=9, font=('Consolas', 9), justify='center')
        self._tr_entry_lo.pack(side=tk.LEFT, padx=(4, 12))
        tk.Label(range_row, text="To", bg=bg, font=('Consolas', 9)).pack(side=tk.LEFT)
        self._tr_entry_hi = tk.Entry(range_row, width=9, font=('Consolas', 9), justify='center')
        self._tr_entry_hi.pack(side=tk.LEFT, padx=(4, 0))
        for entry in (self._tr_entry_lo, self._tr_entry_hi):
            entry.bind('<Return>', self._tr_entry_commit)
            entry.bind('<FocusOut>', self._tr_entry_commit)

        self._tr_range_label = tk.Label(self.f_trim, text="—",
                                        font=('Consolas', 8), bg=bg, fg='#555555')
        self._tr_range_label.pack(pady=(4, 0))

        btn = tk.Frame(self.f_trim, bg=bg)
        btn.pack(fill=tk.X, padx=10, pady=10)
        tk.Button(btn, text="Apply", command=self._tr_applica,
                  bg='#2ecc71', fg='white', font=('Arial', 9, 'bold')).pack(side=tk.LEFT, expand=True, fill=tk.X, padx=2)
        tk.Button(btn, text="Cancel", command=self._tr_annulla).pack(side=tk.LEFT, expand=True, fill=tk.X, padx=2)

    # --- TRIM: apertura pannello ---
    def apri_trim(self):
        sel = self._selezione_effettiva()
        if len(sel) != 1:
            messagebox.showwarning("Trim", "Seleziona esattamente uno spettro.")
            return
        name = self.file_listbox.get(sel[0])
        df = self.spectra[name]['df']
        self._tr_name = name
        self._tr_x = df.index.to_numpy(dtype=float)
        self._tr_y = df[name].to_numpy(dtype=float)
        self._tr_label_nome.config(text=name)

        self.ax.clear()
        self.v_line = self.v_text = None
        tipo = self.spectra[name]['info']['Type']
        if tipo == 'FTIR':
            self.ax.set_xlabel("Wavenumber (cm⁻¹)")
            self.ax.set_ylabel(self._ftir_ylabel([name]))
        elif tipo == 'Raman':
            self.ax.set_xlabel("Raman shift (cm⁻¹)")
            self.ax.set_ylabel("Intensity (A.U.)")
        elif tipo == 'Fluorescence':
            self.ax.set_xlabel("Wavelength (nm)")
            self.ax.set_ylabel("Intensity (A.U.)")
        else:
            self.ax.set_xlabel("Wavelength (nm)")
            self.ax.set_ylabel("Absorbance (A)")
        self.ax.grid(True, linestyle=':', alpha=0.6)

        self._tr_lines['orig'], = self.ax.plot(self._tr_x, self._tr_y,
                                               color='steelblue', alpha=0.4, lw=1.2, label=f'{name} (orig)')
        self._tr_lines['sel'],  = self.ax.plot(self._tr_x, self._tr_y.copy(),
                                               color='green', lw=1.8, label=f'{name}_trim')
        self.ax.legend(fontsize=8)
        if tipo == 'FTIR':
            self.ax.invert_xaxis()

        x_lo, x_hi = float(self._tr_x.min()), float(self._tr_x.max())
        self._tr_vline_lo = self.ax.axvline(x=x_lo, color='purple', ls='--', lw=1.5, alpha=0.8)
        self._tr_vline_hi = self.ax.axvline(x=x_hi, color='purple', ls='--', lw=1.5, alpha=0.8)
        self._tr_drag = None
        self._tr_cids = [
            self.canvas.mpl_connect('button_press_event',   self._tr_drag_start),
            self.canvas.mpl_connect('motion_notify_event',  self._tr_drag_move),
            self.canvas.mpl_connect('button_release_event', self._tr_drag_stop),
        ]
        self._tr_aggiorna_range_label()
        self._tr_aggiorna_preview()
        self.canvas.draw()

        self.f_metadata.pack_forget()
        self.f_trim.pack(fill=tk.BOTH, expand=True)

    # --- TRIM: preview live (evidenzia in verde solo la porzione selezionata) ---
    def _tr_aggiorna_preview(self, *_):
        if self._tr_x is None or 'sel' not in self._tr_lines:
            return
        lo = self._tr_vline_lo.get_xdata()[0]
        hi = self._tr_vline_hi.get_xdata()[0]
        lo, hi = min(lo, hi), max(lo, hi)
        y_sel = np.where((self._tr_x >= lo) & (self._tr_x <= hi), self._tr_y, np.nan)
        self._tr_lines['sel'].set_ydata(y_sel)
        self.canvas.draw_idle()

    def _tr_aggiorna_range_label(self):
        if self._tr_vline_lo is None:
            return
        lo = self._tr_vline_lo.get_xdata()[0]
        hi = self._tr_vline_hi.get_xdata()[0]
        lo, hi = min(lo, hi), max(lo, hi)
        for entry, val in ((self._tr_entry_lo, lo), (self._tr_entry_hi, hi)):
            entry.delete(0, tk.END)
            entry.insert(0, f"{val:.4g}")
        n = int(np.sum((self._tr_x >= lo) & (self._tr_x <= hi)))
        self._tr_range_label.config(text=f"{n} points selected")

    def _tr_snap(self, xv):
        """Aggancia un valore x al punto dati più vicino: i limiti del trim cadono
        sempre su un campione reale dello spettro, mai su un valore interpolato."""
        idx = int(np.argmin(np.abs(self._tr_x - xv)))
        return float(self._tr_x[idx])

    # --- TRIM: inserimento diretto del range da tastiera ---
    def _tr_entry_commit(self, event=None):
        if self._tr_x is None:
            return
        try:
            lo_val = float(self._tr_entry_lo.get())
            hi_val = float(self._tr_entry_hi.get())
        except ValueError:
            self._tr_aggiorna_range_label()   # testo non numerico: ripristina i valori correnti
            return
        lo_val, hi_val = self._tr_snap(min(lo_val, hi_val)), self._tr_snap(max(lo_val, hi_val))
        if lo_val == hi_val:
            idx = int(np.argmin(np.abs(self._tr_x - lo_val)))
            if idx < len(self._tr_x) - 1:
                hi_val = float(self._tr_x[idx + 1])
            elif idx > 0:
                lo_val = float(self._tr_x[idx - 1])
        self._tr_vline_lo.set_xdata([lo_val, lo_val])
        self._tr_vline_hi.set_xdata([hi_val, hi_val])
        self._tr_aggiorna_range_label()
        self._tr_aggiorna_preview()
        self.canvas.draw_idle()

    # --- TRIM: drag delle linee verticali (agganciate al passo reale dei dati) ---
    def _tr_drag_start(self, event):
        if event.inaxes != self.ax or event.button != 1:
            return
        for key, line in [('lo', self._tr_vline_lo), ('hi', self._tr_vline_hi)]:
            x_disp = self.ax.transData.transform((line.get_xdata()[0], 0))[0]
            if abs(event.x - x_disp) < 8:
                self._tr_drag = key
                return

    def _tr_drag_move(self, event):
        if self._tr_drag is None or event.inaxes != self.ax or event.xdata is None:
            return
        xv = self._tr_snap(event.xdata)
        lo = self._tr_vline_lo.get_xdata()[0]
        hi = self._tr_vline_hi.get_xdata()[0]
        if self._tr_drag == 'lo':
            if xv >= hi:
                return
            self._tr_vline_lo.set_xdata([xv, xv])
        else:
            if xv <= lo:
                return
            self._tr_vline_hi.set_xdata([xv, xv])
        self._tr_aggiorna_range_label()
        self._tr_aggiorna_preview()
        self.canvas.draw_idle()

    def _tr_drag_stop(self, event):
        self._tr_drag = None

    # --- TRIM: applica (crea una nuova traccia con la sola porzione selezionata) ---
    def _tr_applica(self):
        lo = self._tr_vline_lo.get_xdata()[0]
        hi = self._tr_vline_hi.get_xdata()[0]
        lo, hi = min(lo, hi), max(lo, hi)
        mask = (self._tr_x >= lo) & (self._tr_x <= hi)
        if mask.sum() < 2:
            messagebox.showwarning("Trim", "La selezione contiene meno di 2 punti.")
            return
        new_name = f"{self._tr_name}_trim"
        final_name = new_name
        counter = 2
        while final_name in self.spectra:
            final_name = f"{new_name}_{counter}"
            counter += 1
        new_info = dict(self.spectra[self._tr_name]['info'])
        new_info['Sample']       = final_name
        new_info['Derived from'] = self._tr_name
        new_info['Operation']    = 'Trim'
        new_info['Trim Range']   = f"{lo:.2f} – {hi:.2f}"
        orig_index = self.spectra[self._tr_name]['df'].index
        self.spectra[final_name] = {
            'df':   pd.DataFrame({final_name: self._tr_y[mask]}, index=orig_index[mask]),
            'info': new_info,
        }
        self._dirty = True
        self._tr_chiudi_pannello()
        self.aggiorna_vista()

    # --- TRIM: annulla senza modifiche ---
    def _tr_annulla(self):
        self._tr_chiudi_pannello()
        self.aggiorna_vista()

    # --- TRIM: chiude il pannello e ripristina metadata ---
    def _tr_chiudi_pannello(self):
        for cid in self._tr_cids:
            self.canvas.mpl_disconnect(cid)
        self._tr_cids     = []
        self._tr_drag     = None
        self._tr_vline_lo = None
        self._tr_vline_hi = None
        self._tr_name     = None
        self._tr_x        = None
        self._tr_y        = None
        self._tr_lines    = {}
        self.f_trim.pack_forget()
        self.f_metadata.pack(fill=tk.BOTH, expand=True)

    # --- CURSORE ---
    def _cursor_panel_attivo(self):
        return (self._sc_name is not None or self._sub_name_a is not None or self._norm_names
                or self._bl_name is not None or self._ab_name is not None
                or self._sm_name is not None or self._dc_name is not None
                or self._tr_name is not None)

    def _cursor_best_position(self):
        """Sceglie l'angolo del grafico con minor densità di punti tracciati nelle
        vicinanze, per evitare che il riquadro X/Y copra i dati (stile 'best' delle
        legende di matplotlib). Restituisce (ax_x, ax_y, ha, va) in coordinate assi."""
        corners = [(0.02, 0.98, 'left', 'top'), (0.98, 0.98, 'right', 'top'),
                   (0.02, 0.02, 'left', 'bottom'), (0.98, 0.02, 'right', 'bottom')]
        xlim, ylim = self.ax.get_xlim(), self.ax.get_ylim()
        dx, dy = xlim[1] - xlim[0], ylim[1] - ylim[0]
        if dx <= 0 or dy <= 0 or not self._cursor_data:
            return corners[0]
        best, best_count = corners[0], None
        for cx, cy, ha, va in corners:
            x0, x1 = (xlim[0], xlim[0] + 0.35 * dx) if cx < 0.5 else (xlim[1] - 0.35 * dx, xlim[1])
            y0, y1 = (ylim[0], ylim[0] + 0.35 * dy) if cy < 0.5 else (ylim[1] - 0.35 * dy, ylim[1])
            count = sum(int(np.count_nonzero((xs >= x0) & (xs <= x1) & (ys >= y0) & (ys <= y1)))
                        for _, xs, ys in self._cursor_data)
            if best_count is None or count < best_count:
                best, best_count = (cx, cy, ha, va), count
        return best

    def _cursor_drag_start(self, event):
        if self.v_text is None or event.button != 1 or self._cursor_panel_attivo():
            return
        renderer = self.canvas.get_renderer()
        if renderer is None:
            return
        if self.v_text.get_window_extent(renderer).contains(event.x, event.y):
            self._cursor_drag = True

    def _cursor_drag_stop(self, event):
        self._cursor_drag = False

    def on_mouse_move(self, event):
        if self._cursor_drag:
            if self.v_text is not None and event.x is not None and event.y is not None:
                ax_x, ax_y = self.ax.transAxes.inverted().transform((event.x, event.y))
                ha = 'left' if ax_x <= 0.5 else 'right'
                va = 'bottom' if ax_y <= 0.5 else 'top'
                self._cursor_box_pos = (ax_x, ax_y, ha, va)
                self.v_text.set_position((ax_x, ax_y))
                self.v_text.set_ha(ha)
                self.v_text.set_va(va)
                self.canvas.draw_idle()
            return
        if not event.inaxes or not self.spectra or self._cursor_panel_attivo():
            return
        x = event.xdata
        if self.v_line: self.v_line.remove()
        if self.v_text: self.v_text.remove()
        self.v_line = self.ax.axvline(x=x, color='red', linestyle='--', lw=0.7)
        cursor_txt = f"X: {x:.2f}\n"
        for name, xs, ys in self._cursor_data:
            idx = np.argmin(np.abs(xs - x))
            cursor_txt += f"{name[:15]}: {ys[idx]:.4f}\n"
        ax_x, ax_y, ha, va = self._cursor_box_pos or self._cursor_best_position()
        self.v_text = self.ax.text(
            ax_x, ax_y, cursor_txt, transform=self.ax.transAxes,
            horizontalalignment=ha, verticalalignment=va, fontsize=8,
            bbox=dict(facecolor='white', alpha=0.7)
        )
        self.canvas.draw_idle()

    def on_mouse_leave(self, event):
        if self._cursor_drag:
            return
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

    def _menu_grafico(self, event):
        """Menu contestuale sul grafico (tasto destro)."""
        if self._dc_name is not None:
            return   # in deconvoluzione il tasto destro rimuove i picchi
        menu = tk.Menu(self.root, tearoff=0)
        stato = tk.DISABLED if self._pannello_attivo() else tk.NORMAL
        menu.add_command(label="Paste from clipboard   (Ctrl+V)",
                         command=self._incolla_clipboard, state=stato)
        menu.tk_popup(event.x_root, event.y_root)

    def _deduci_tipo_incolla(self, header_x):
        """Deduce il tipo dall'intestazione X; None se non chiaro (→ chiedi).
        Nota: il numero d'onda (cm⁻¹) è usato sia da FTIR che da Raman, quindi
        da solo NON è discriminante → si chiede all'utente."""
        h = (header_x or '').lower()
        if 'raman' in h:
            return 'Raman'
        if 'ftir' in h:
            return 'FTIR'
        if 'wavelength' in h or 'nm' in h:
            return 'UV-Vis'
        return None

    def _chiedi_tipo_spettro(self):
        """Dialogo modale: chiede il tipo di spettro. Restituisce il tipo o None."""
        win = tk.Toplevel(self.root)
        win.title("Tipo di spettro")
        win.transient(self.root)
        win.resizable(False, False)
        tk.Label(win, text="Tipo dello spettro incollato:",
                 font=('Arial', 10)).pack(padx=20, pady=(15, 10))
        scelta = {'val': None}

        def pick(t):
            scelta['val'] = t
            win.destroy()

        fr = tk.Frame(win)
        fr.pack(padx=20, pady=(0, 15))
        for t in ('UV-Vis', 'FTIR', 'Raman'):
            tk.Button(fr, text=t, width=9, command=lambda tt=t: pick(tt)).pack(side=tk.LEFT, padx=4)
        win.update_idletasks()
        x = self.root.winfo_rootx() + (self.root.winfo_width() - win.winfo_width()) // 2
        y = self.root.winfo_rooty() + (self.root.winfo_height() - win.winfo_height()) // 2
        win.geometry(f"+{x}+{y}")
        win.grab_set()
        self.root.wait_window(win)
        return scelta['val']

    def _incolla_clipboard(self, event=None):
        """Incolla due colonne x,y dalla clipboard come nuovo spettro (Ctrl+V sul grafico)."""
        if self._pannello_attivo():
            return
        try:
            text = self.root.clipboard_get()
        except tk.TclError:
            return
        lines = [ln for ln in text.splitlines() if ln.strip()]
        if not lines:
            return

        # Rileva separatore: tab (Origin/Excel) → ; (europeo, decimale virgola) → spazi → ,
        sample = lines[0]
        if '\t' in text:
            sep, dec_comma = '\t', False
        elif ';' in sample:
            sep, dec_comma = ';', True
        elif re.search(r'\S\s+\S', sample):
            sep, dec_comma = None, False   # split su spazi
        else:
            sep, dec_comma = ',', False

        def to_float(s):
            return float(s.replace(',', '.') if dec_comma else s)

        # Intestazione: se la prima riga non è numerica, usala per nome ed etichetta X
        header_x = header_name = None
        fparts = [p for p in (lines[0].split(sep) if sep else lines[0].split()) if p != '']
        if len(fparts) >= 2:
            try:
                to_float(fparts[0]); to_float(fparts[1])
            except ValueError:
                header_x, header_name = fparts[0].strip(), fparts[1].strip()

        xs, ys = [], []
        for ln in lines:
            parts = ln.split(sep) if sep else ln.split()
            parts = [p for p in parts if p != '']
            if len(parts) < 2:
                continue
            try:
                xv, yv = to_float(parts[0]), to_float(parts[1])
            except ValueError:
                continue   # riga di intestazione o non numerica → saltata
            xs.append(xv)
            ys.append(yv)

        if len(xs) < 2:
            messagebox.showwarning("Paste", "La clipboard non contiene due colonne numeriche valide.")
            return

        x = np.array(xs)
        y = self._round_sig(np.array(ys))
        # Tipo: dall'intestazione se chiaro, altrimenti chiedi all'utente
        tipo = self._deduci_tipo_incolla(header_x)
        if tipo is None:
            tipo = self._chiedi_tipo_spettro()
            if tipo is None:
                return   # annullato
        base = header_name if header_name else 'pasted'
        final = base
        counter = 2
        while final in self.spectra:
            final = f"{base}_{counter}"
            counter += 1
        info = {'Sample': final, 'Date': 'N/A', 'Type': tipo, 'Source': 'Clipboard'}
        if header_x:
            info['X Axis'] = header_x
        self.spectra[final] = {
            'df':   pd.DataFrame({final: y}, index=x),
            'info': info,
        }
        self._dirty = True   # incollato dalla clipboard: non esiste su disco, perderlo sarebbe definitivo
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
            traceback.print_exc()
            messagebox.showerror("Error", f"Could not read {os.path.basename(path)}: {e}")

    @staticmethod
    def _round_sig(arr, sig=4):
        """Arrotonda un array numpy a `sig` cifre significative (vettorizzato)."""
        arr = np.asarray(arr, dtype=float)
        out = arr.copy()
        mask = np.isfinite(arr) & (arr != 0)
        vals = arr[mask]
        scale = 10.0 ** (sig - 1 - np.floor(np.log10(np.abs(vals))))
        out[mask] = np.round(vals * scale) / scale
        return out

    def leggi_dsp(self, path):
        with open(path, 'r', encoding='latin-1') as f:
            lines = [line.strip() for line in f.readlines()]
        name = os.path.basename(path).replace('.dsp', '')
        wl_start, wl_end, n_punti = float(lines[5]), float(lines[6]), int(lines[8])
        idx_data = lines.index("#DATA") + 1
        y_vals = self._round_sig([float(lines[i]) for i in range(idx_data, idx_data + n_punti)])
        x_vals = np.linspace(wl_start, wl_end, n_punti)
        self.spectra[name] = {
            'df': pd.DataFrame({'x': x_vals, name: y_vals}).set_index('x'),
            'info': {'Sample': name, 'Date': lines[16], 'Type': 'UV-Vis'}
        }

    def leggi_sp(self, path):
        """Parser nativo per file binari .SP di PerkinElmer FTIR.

        Struttura file: magic "PEPE" (4B) + descrizione (40B) + blocchi top-level.
        Il blocco id=120 (DSet2DC1DI) contiene member-block con i dati spettrali.
        Member id rilevanti (int16 little-endian, negativi):
          -29838  AbscissaRange        : 2 double → StartX, EndX
          -29835  NumPoints            : Int32    → numero punti
          -29833  XAxisLabel           : stringa  → unità asse X
          -29832  YAxisLabel           : stringa  → tipo misura (A/T/R)
          -29828  DataSetData          : Int32 (byte) + n double Y
          -29827  DataSetName          : stringa  → percorso file originale
          -29825  DataSetHistoryRecord : testo    → data, operatore, strumento, …
        """
        with open(path, 'rb') as f:
            raw = f.read()

        if raw[:4] != b'PEPE':
            raise ValueError("File .SP non valido: magic 'PEPE' mancante")

        DSET2DC1DI     = 120
        ABSCISSA_RANGE = -29838
        NUM_POINTS     = -29835
        X_LABEL        = -29833
        Y_LABEL        = -29832
        DATA           = -29828
        ORIG_NAME      = -29827
        HISTORY        = -29825

        x_start = x_end = n_points = None
        y_data    = None
        x_label   = None
        y_label   = None
        orig_name = None
        meta      = {}

        pos = 44  # salta magic (4B) + descrizione (40B)
        while pos + 6 <= len(raw):
            block_id  = struct.unpack_from('<h', raw, pos)[0]
            block_len = struct.unpack_from('<i', raw, pos + 2)[0]
            pos += 6
            block_data = raw[pos:pos + block_len]
            pos += block_len

            if block_id != DSET2DC1DI:
                continue

            mpos = 0
            while mpos + 8 <= len(block_data):
                m_id  = struct.unpack_from('<h', block_data, mpos)[0]
                m_len = struct.unpack_from('<i', block_data, mpos + 2)[0]
                mpos += 6
                m_data = block_data[mpos + 2:mpos + m_len]  # salta TypeCode (2B)
                mpos  += m_len

                if m_id == HISTORY:
                    meta = self._sp_extract_metadata(m_data)
                elif m_id == ABSCISSA_RANGE and len(m_data) >= 16:
                    x_start = struct.unpack_from('<d', m_data, 0)[0]
                    x_end   = struct.unpack_from('<d', m_data, 8)[0]
                elif m_id == NUM_POINTS and len(m_data) >= 4:
                    n_points = struct.unpack_from('<i', m_data, 0)[0]
                elif m_id == DATA and len(m_data) >= 4:
                    byte_len = struct.unpack_from('<i', m_data, 0)[0]
                    n_elem   = byte_len // 8
                    y_data   = np.array(struct.unpack_from(f'<{n_elem}d', m_data, 4))
                elif m_id == X_LABEL:
                    x_label = self._sp_read_string(m_data)
                elif m_id == Y_LABEL:
                    y_label = self._sp_read_string(m_data)
                elif m_id == ORIG_NAME:
                    orig_name = self._sp_read_string(m_data)
            break

        if x_start is None or x_end is None or y_data is None:
            raise ValueError("Struttura .SP non riconosciuta: campi dati mancanti")

        if n_points is None:
            n_points = len(y_data)

        Y_LABEL_MAP = {
            'A': 'Absorbance', 'T': 'Transmittance', 'R': 'Reflectance',
            '%T': 'Transmittance %', '%R': 'Reflectance %', 'ATR': 'ATR',
        }

        x_vals = np.linspace(x_start, x_end, n_points)
        name   = os.path.splitext(os.path.basename(path))[0]
        info = {
            'Sample':     name,
            'Date':       meta.get('Date', 'N/A'),
            'Type':       'FTIR',
            'Y Axis':     Y_LABEL_MAP.get(y_label, y_label or 'N/A'),
            'X Axis':     x_label or 'cm-1',
        }
        for key in ('Operator', 'Instrument', 'Detector', 'Apodization', 'Firmware', 'Processing'):
            if key in meta:
                info[key] = meta[key]
        if orig_name:
            info['Original File'] = os.path.basename(orig_name)

        self.spectra[name] = {
            'df':   pd.DataFrame({'x': x_vals, name: self._round_sig(y_data)}).set_index('x'),
            'info': info
        }

    def _ftir_ylabel(self, names=None):
        """Restituisce l'etichetta asse Y per spettri FTIR basandosi sui metadati.

        Se tutti gli spettri FTIR selezionati hanno la stessa unità Y la usa;
        altrimenti cade su 'Transmittance / Absorbance'.
        """
        UNITS = {
            'Absorbance':     'Absorbance (A)',
            'Transmittance':  'Transmittance (T)',
            'Transmittance %':'Transmittance (%T)',
            'Reflectance':    'Reflectance (R)',
            'Reflectance %':  'Reflectance (%R)',
            'ATR':            'ATR Absorbance',
        }
        if names is None:
            names = list(self.spectra.keys())
        labels = {
            self.spectra[n]['info'].get('Y Axis', '')
            for n in names
            if self.spectra.get(n, {}).get('info', {}).get('Type') == 'FTIR'
        } - {'', 'N/A'}
        if len(labels) == 1:
            return UNITS.get(labels.pop(), 'Transmittance / Absorbance')
        return 'Transmittance / Absorbance'

    @staticmethod
    def _sp_read_string(data):
        """Legge una stringa in formato PerkinElmer: Int16 length + bytes ASCII."""
        if len(data) < 2:
            return None
        slen = struct.unpack_from('<h', data, 0)[0]
        if slen <= 0 or len(data) < 2 + slen:
            return None
        return data[2:2 + slen].decode('ascii', errors='replace').strip()

    @staticmethod
    def _sp_extract_metadata(history_data):
        """Estrae metadati dal DataSetHistoryRecord.

        Il blocco contiene una sequenza di record (operatore, operazione, timestamp,
        parametri) separati da byte-tipo (0x73/0x74/0x75/0x76). Ogni stringa leggibile
        è in formato latin-1 e i campi utili sono estratti con regex mirate.
        """
        meta = {}
        try:
            text = history_data.decode('latin-1', errors='replace')

            # Data originale: voce "Created by Instrument" + timestamp ctime
            m = re.search(
                r'Created by Instrument.{1,10}([A-Za-z]{3} [A-Za-z]{3} {1,2}\d{1,2} \d{2}:\d{2}:\d{2} \d{4})',
                text, re.DOTALL
            )
            if m:
                try:
                    dt = datetime.datetime.strptime(m.group(1).strip(), '%a %b %d %H:%M:%S %Y')
                    meta['Date'] = dt.strftime('%Y-%m-%d %H:%M')
                except Exception:
                    pass

            # Operatore: stringa leggibile subito prima di "Created by Instrument"
            ci_idx = text.find('Created by Instrument')
            if ci_idx > 0:
                before = text[max(0, ci_idx - 60):ci_idx]
                words = re.findall(r'[A-Za-z][A-Za-z0-9 _.-]{1,28}', before)
                if words:
                    meta['Operator'] = words[-1].strip()

            # Modello strumento (es. "Spectrum Series FTIR", "Spectrum Two")
            m = re.search(
                r'(Spectrum(?:\s+Series)?\s+(?:FTIR|UV|GX|One|Two|100|65|400|BX)[\w\s]{0,20})',
                text
            )
            if m:
                meta['Instrument'] = m.group(1).strip()

            # Firmware / versione software
            m = re.search(r'(CPU\d+\s+\w+\s+[\d.]+\s+[\d\w-]+\s+[\d:]+)', text)
            if m:
                meta['Firmware'] = m.group(1).strip()

            # Detector
            for det in ['MIR TGS', 'MIR MCT', 'NIR InGaAs', 'MCT/A', 'DTGS', 'MCT', 'TGS', 'InGaAs']:
                if det in text:
                    meta['Detector'] = det
                    break

            # Apodizzazione
            for apod in ['Norton-Beer Strong', 'Norton-Beer Medium', 'Norton-Beer Weak',
                         'Happ-Genzel', 'Blackman-Harris', 'Strong', 'Medium', 'Weak',
                         'Cosine', 'Triangular', 'Boxcar']:
                if apod in text:
                    meta['Apodization'] = apod
                    break

            # Storia elaborazioni (ordine di comparsa, senza duplicati)
            ops = re.findall(
                r'(Absorbance|Transmittance|Reflectance|ATR Correction|'
                r'Atmospheric Correction|Baseline Correction|Normalize|Smooth|'
                r'Add|Subtract|Div|Multiply)',
                text
            )
            seen = set()
            unique = [x for x in ops if not (x in seen or seen.add(x))]
            if unique:
                meta['Processing'] = ' → '.join(unique)

        except Exception:
            pass
        return meta

    # Pattern per riconoscere intestazioni che sono sole unità di misura
    _UNIT_RE = re.compile(
        r'^(nm|cm\s*[-–^]?\s*1|a\.?u\.?|abs(orbance)?|transmittance|'
        r'intensity|int|wavelength|wavenumber|[atr])$',
        re.IGNORECASE
    )

    def leggi_csv(self, path):
        with open(path, 'r', encoding='latin-1') as f:
            head = f.read(4096)
        first_lines = head.splitlines()

        # Report Cary di soli parametri (nessun dato spettrale)
        if 'Scan Analysis Report' in head:
            messagebox.showinfo(
                "Cary Report",
                f"'{os.path.basename(path)}' è un report di soli parametri Cary "
                "(Scan Analysis Report): contiene i metadati strumentali ma nessuno "
                "spettro. I dati numerici vanno esportati con 'File Storage' attivo.")
            return

        line1 = first_lines[1] if len(first_lines) > 1 else ''
        if 'Wavelength' in line1 and ('Abs' in line1 or 'abs' in line1):
            self._leggi_csv_cary(path)
            return
        df = pd.read_csv(path, index_col=0)
        tipo = 'FTIR' if df.index.max() > 2000 else 'UV-Vis'
        filename_base = os.path.splitext(os.path.basename(path))[0]

        # Se tutte le colonne hanno nomi che sono solo unità di misura,
        # usa il nome del file come etichetta base
        all_units = all(self._UNIT_RE.match(str(c).strip()) for c in df.columns)

        for col in df.columns:
            serie = df[[col]].dropna()

            if all_units:
                label = filename_base if len(df.columns) == 1 else f"{filename_base}_{col}"
            else:
                label = str(col)

            # Evita collisioni con spettri già caricati
            final_label = label
            counter = 2
            while final_label in self.spectra:
                final_label = f"{label}_{counter}"
                counter += 1

            # Il nome della colonna nel DataFrame deve coincidere con la chiave
            serie = serie.rename(columns={col: final_label})
            serie[final_label] = self._round_sig(serie[final_label].to_numpy())

            self.spectra[final_label] = {
                'df': serie,
                'info': {'Sample': final_label, 'Date': 'N/A', 'Type': tipo}
            }

    def _leggi_csv_cary(self, path):
        with open(path, 'r', encoding='latin-1') as f:
            lines = f.readlines()

        # Rileva il separatore dalla riga di intestazione
        sep = ';' if ';' in lines[1] else ','

        # Riga 0: nomi campioni — ogni spettro occupa 2 colonne, quindi i nomi
        # sono separati da doppio separatore; le virgole nel nome sopravvivono
        raw_names = [s.strip() for s in lines[0].split(sep * 2) if s.strip()]
        filename_base = os.path.splitext(os.path.basename(path))[0]

        # Riga 1: intestazioni — conta quante coppie (Wavelength, Abs) ci sono
        headers = [s.strip() for s in lines[1].split(sep)]
        n_spectra = sum(1 for h in headers if 'wavelength' in h.lower())
        if n_spectra == 0:
            n_spectra = 1

        # Parsing riga per riga — accesso posizionale alle colonne, così spettri
        # di lunghezza diversa (celle vuote intermedie) non slittano di posto
        columns = [[] for _ in range(n_spectra * 2)]
        for line in lines[2:]:
            parts = [s.strip() for s in line.split(sep)]
            for i in range(n_spectra):
                xi, yi = i * 2, i * 2 + 1
                if yi < len(parts) and parts[xi] and parts[yi]:
                    try:
                        columns[xi].append(float(parts[xi]))
                        columns[yi].append(float(parts[yi]))
                    except ValueError:
                        pass

        # Estrai metadati dal blocco in fondo al file
        meta = {}
        for line in lines[2:]:
            stripped = line.strip()
            for key, pattern in [('Date',       'Collection Time:'),
                                  ('Instrument', 'Instrument '),
                                  ('Start (nm)', 'Start (nm)'),
                                  ('Stop (nm)',  'Stop (nm)'),
                                  ('Interval',   'UV-Vis Data Interval'),
                                  ('SBW (nm)',   'UV-Vis SBW')]:
                if stripped.startswith(pattern) and key not in meta:
                    val = stripped[len(pattern):].split(sep)[0].strip().lstrip(':').strip()
                    if val:
                        meta[key] = val

        for i in range(n_spectra):
            x = np.array(columns[i * 2])
            y = self._round_sig(np.array(columns[i * 2 + 1]))
            if len(x) == 0:
                continue

            label = raw_names[i] if i < len(raw_names) and raw_names[i] \
                    else (filename_base if n_spectra == 1 else f"{filename_base}_{i+1}")
            final_label = label
            counter = 2
            while final_label in self.spectra:
                final_label = f"{label}_{counter}"
                counter += 1

            tipo = 'FTIR' if x.max() > 2000 else 'UV-Vis'
            info = {'Sample': final_label, 'Date': meta.get('Date', 'N/A'), 'Type': tipo}
            info.update({k: v for k, v in meta.items() if k != 'Date'})
            self.spectra[final_label] = {
                'df':   pd.DataFrame({final_label: y}, index=x),
                'info': info,
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
        y_data = self._round_sig(np.frombuffer(y_segment[:n_points * 4], dtype=np.float32).copy())
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
                'File':            os.path.basename(path),
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
            self._cursor_data = []
            self.canvas.draw()
            return

        for idx, sname in enumerate(self.spectra.keys()):
            self.file_listbox.insert(tk.END, sname)
            if sname in self._hidden_spectra:
                self.file_listbox.itemconfig(idx, fg='#aaaaaa')

        # Con un solo spettro caricato, selezionalo subito: evita di dover cliccare
        # sull'unica voce prima di usare i pulsanti che richiedono una selezione.
        if len(self.spectra) == 1:
            self.file_listbox.selection_set(0)

        # Cache array numpy per il cursore (evita overhead pandas a ogni mouse move);
        # gli spettri nascosti non compaiono nel grafico quindi neanche nel cursore.
        self._cursor_data = [
            (name, d['df'].index.to_numpy(dtype=float), d['df'][name].to_numpy(dtype=float))
            for name, d in self.spectra.items() if name not in self._hidden_spectra
        ]

        master_df = pd.concat([s['df'] for s in self.spectra.values()], axis=1).sort_index()

        # Plotting — rilevamento tipo dal campo 'Type'
        types = {d['info']['Type'] for d in self.spectra.values()}

        # Tabella tabulata — etichetta asse X coerente col tipo di spettro
        if 'FTIR' in types:
            x_index_label = 'Wavenumber (cm-1)'
        elif 'Raman' in types:
            x_index_label = 'Raman shift (cm-1)'
        else:
            x_index_label = 'Wavelength (nm)'
        table_df = master_df.copy()
        table_df.index.name = x_index_label
        self.data_box.insert(tk.END, table_df.round(4).to_csv(sep='\t', lineterminator='\n'))

        visible_cols = [c for c in master_df.columns if c not in self._hidden_spectra]
        for col in visible_cols:
            valid = master_df[col].dropna()
            self.ax.plot(valid.index, valid.values, label=col, lw=1.5)

        if visible_cols:
            self.ax.set_xlim(min(self.spectra[c]['df'].index.min() for c in visible_cols),
                             max(self.spectra[c]['df'].index.max() for c in visible_cols))
        else:
            self.ax.set_xlim(master_df.index.min(), master_df.index.max())

        if 'FTIR' in types:
            self.ax.set_xlabel("Wavenumber (cm⁻¹)")
            self.ax.set_ylabel(self._ftir_ylabel())
            self.ax.invert_xaxis()
        elif 'Raman' in types:
            self.ax.set_xlabel("Raman shift (cm⁻¹)")
            self.ax.set_ylabel("Intensity (A.U.)")
        elif 'Fluorescence' in types:
            self.ax.set_xlabel("Wavelength (nm)")
            self.ax.set_ylabel("Intensity (A.U.)")
        else:
            self.ax.set_xlabel("Wavelength (nm)")
            self.ax.set_ylabel("Absorbance (A)")

        if visible_cols:
            self.ax.legend(fontsize='8', loc='best')
        self.ax.grid(True, linestyle=':', alpha=0.6)
        self.canvas.draw()

        # Metadata panel — display generico di tutte le chiavi info
        for name, data in self.spectra.items():
            self.param_box.insert(tk.END, f"{'='*38}\n")
            for k, v in data['info'].items():
                self.param_box.insert(tk.END, f"{k:<20}: {v}\n")
            self.param_box.insert(tk.END, "\n")

    # --- OPERAZIONI ARITMETICHE ---
    def media_selezionati(self):
        sel = self.file_listbox.curselection()
        if len(sel) < 2:
            messagebox.showwarning("Average", "Seleziona almeno 2 spettri.")
            return
        names = [self.file_listbox.get(i) for i in sel]

        # Verifica che i tipi siano compatibili
        types = {self.spectra[n]['info']['Type'] for n in names}
        if len(types) > 1:
            messagebox.showwarning("Average", f"Gli spettri selezionati hanno tipi diversi: {', '.join(types)}.")
            return

        # Calcola la media interpolando sull'indice comune (unione)
        dfs = [self.spectra[n]['df'].rename(columns={n: n}) for n in names]
        combined = pd.concat(dfs, axis=1)
        mean_series = combined.mean(axis=1)

        # Nome del risultato
        short_names = '_'.join(n[:8] for n in names)
        avg_label = f"avg({short_names})"
        # Evita collisioni
        final_label = avg_label
        counter = 2
        while final_label in self.spectra:
            final_label = f"{avg_label}_{counter}"
            counter += 1

        avg_df = pd.DataFrame({final_label: mean_series})
        avg_df.index.name = combined.index.name

        self.spectra[final_label] = {
            'df': avg_df,
            'info': {
                'Sample':      final_label,
                'Date':        'N/A',
                'Type':        next(iter(types)),
                'Operation':   'Average',
                'Source':      ', '.join(names),
            }
        }
        self._dirty = True
        self.aggiorna_vista()

    # --- CONTEXT MENU LISTBOX ---
    def _listbox_tasto_destro(self, event):
        idx = self.file_listbox.nearest(event.y)
        if idx < 0 or idx >= self.file_listbox.size():
            return
        name = self.file_listbox.get(idx)
        # Se l'elemento cliccato fa parte di una selezione multipla, opera su tutta
        sel = self.file_listbox.curselection()
        names = [self.file_listbox.get(i) for i in sel] if idx in sel else [name]

        menu = tk.Menu(self.root, tearoff=0)
        menu.add_command(label="Copy X,Y values to clipboard",
                         command=lambda: self._copia_xy_clipboard(name))
        etichetta = "Copy for Origin (Long Name/Units/Comment)" if len(names) == 1 \
            else f"Copy {len(names)} spectra for Origin"
        menu.add_command(label=etichetta,
                         command=lambda: self._copia_origin(names))
        menu.add_separator()
        menu.add_command(label="Hide" if len(names) == 1 else f"Hide {len(names)} spectra",
                         command=lambda: self._nascondi_spettri(names))
        menu.add_command(label="Show" if len(names) == 1 else f"Show {len(names)} spectra",
                         command=lambda: self._mostra_spettri(names))
        menu.tk_popup(event.x_root, event.y_root)

    def _nascondi_spettri(self, names):
        self._hidden_spectra.update(names)
        self.aggiorna_vista()

    def _mostra_spettri(self, names):
        self._hidden_spectra.difference_update(names)
        self.aggiorna_vista()

    def _copia_xy_clipboard(self, name):
        df = self.spectra[name]['df']
        testo = "\n".join(f"{x}\t{y}" for x, y in zip(df.index, df[name]))
        self.root.clipboard_clear()
        self.root.clipboard_append(testo)

    def _origin_headers(self, name):
        """Restituisce (x_long, x_unit, y_long, y_unit) secondo il tipo di spettro."""
        info = self.spectra[name]['info']
        tipo = info.get('Type', 'UV-Vis')
        if tipo == 'FTIR':
            y_long = info.get('Y Axis', 'Absorbance')
            y_unit = {'Transmittance': '%T', 'Transmittance %': '%T',
                      'Reflectance': '%R', 'Reflectance %': '%R'}.get(y_long, '')
            return 'Wavenumber', 'cm-1', y_long, y_unit
        if tipo == 'Raman':
            return 'Raman shift', 'cm-1', 'Intensity', 'A.U.'
        if tipo == 'Fluorescence':
            return 'Wavelength', 'nm', 'Intensity', 'A.U.'
        return 'Wavelength', 'nm', 'Absorbance', ''

    def _copia_origin(self, names):
        if isinstance(names, str):
            names = [names]

        # Allinea per indice: se la X è identica risulta una sola colonna condivisa,
        # altrimenti l'unione degli indici con celle vuote dove un dato manca
        merged = pd.concat([self.spectra[n]['df'] for n in names], axis=1).sort_index()

        # Intestazioni: X dal primo spettro, poi una Y per ogni spettro
        x_long, x_unit, _, _ = self._origin_headers(names[0])
        long_row = [x_long]
        unit_row = [x_unit]
        comm_row = ['']
        for n in names:
            _, _, y_long, y_unit = self._origin_headers(n)
            info = self.spectra[n]['info']
            # FTIR: usa 'Sample' (il 'Original File' interno è troncato in 8.3 DOS)
            if info.get('Type') == 'FTIR':
                comment = info.get('Sample') or n
            else:
                comment = info.get('File') or info.get('Original File') or n
            long_row.append(y_long)
            unit_row.append(y_unit)
            comm_row.append(comment)

        righe = ['\t'.join(long_row), '\t'.join(unit_row), '\t'.join(comm_row)]
        for idx, row in zip(merged.index, merged[names].to_numpy()):
            celle = [f"{v}" if pd.notna(v) else '' for v in row]
            righe.append(f"{idx}\t" + '\t'.join(celle))

        self.root.clipboard_clear()
        self.root.clipboard_append("\n".join(righe))

    # --- GESTIONE LISTA ---
    def _selezione_effettiva(self):
        """Selezione corrente; se nulla ma c'è un solo spettro, lo auto-seleziona."""
        sel = self.file_listbox.curselection()
        if not sel and self.file_listbox.size() == 1:
            self.file_listbox.selection_set(0)
            return (0,)
        return sel

    def _pannello_attivo(self):
        """Restituisce il nome del pannello operativo aperto, o None."""
        if self._sc_name is not None:
            return "Scattering Correction"
        if self._sub_name_a is not None:
            return "Scaled Subtraction"
        if self._norm_names:
            return "Normalize"
        if self._bl_name is not None:
            return "Baseline (FTIR)"
        if self._ab_name is not None:
            return "Adaptive Baseline"
        if self._sm_name is not None:
            return "Smoothing"
        if self._dc_name is not None:
            return "Deconvolution"
        if self._tr_name is not None:
            return "Trim"
        return None

    def remove_selected(self):
        pannello = self._pannello_attivo()
        if pannello:
            messagebox.showwarning("Remove", f"Chiudi prima il pannello {pannello} (Apply o Cancel).")
            return
        selected = self.file_listbox.curselection()
        for i in reversed(selected):
            name = self.file_listbox.get(i)
            if name in self.spectra: del self.spectra[name]
            self._hidden_spectra.discard(name)
        self.aggiorna_vista()

    def clear_all(self):
        pannello = self._pannello_attivo()
        if pannello:
            messagebox.showwarning("Clear All", f"Chiudi prima il pannello {pannello} (Apply o Cancel).")
            return
        self.spectra = {}
        self._hidden_spectra = set()
        self._dirty = False
        self.aggiorna_vista()

    # --- MOSTRA/NASCONDI: escludono uno spettro dal grafico senza rimuoverlo ---
    def hide_selected(self):
        sel = self.file_listbox.curselection()
        if not sel:
            return
        for i in sel:
            self._hidden_spectra.add(self.file_listbox.get(i))
        self.aggiorna_vista()

    def show_selected(self):
        sel = self.file_listbox.curselection()
        if not sel:
            return
        for i in sel:
            self._hidden_spectra.discard(self.file_listbox.get(i))
        self.aggiorna_vista()

    def esporta_csv(self):
        if not self.spectra: return
        path = filedialog.asksaveasfilename(defaultextension=".csv", initialdir=self.current_dir)
        if path:
            dfs = [s['df'] for s in self.spectra.values()]
            pd.concat(dfs, axis=1).round(4).to_csv(path)
            self._dirty = False

    def on_exit(self):
        """Chiude l'app, avvisando se ci sono risultati calcolati (scattering, fit,
        medie, trim, ...) non ancora esportati con 'Export CSV'."""
        if self._dirty and not messagebox.askyesno(
                "Exit",
                "You have unsaved calculated results (not exported via Export CSV). Exit anyway?"):
            return
        self.root.destroy()

    def _pulisci_cursore(self):
        """Rimuove l'overlay interattivo del cursore (linea verticale + testo) dagli
        assi, così non finisce incollato nella figura salvata/passata all'editor."""
        if self.v_line is not None:
            self.v_line.remove()
            self.v_line = None
        if self.v_text is not None:
            self.v_text.remove()
            self.v_text = None

    def salva_figura_pickle(self):
        """Salva la figura matplotlib corrente come pickle: riapribile con
        plot_editor.pyw come oggetto Figure/Axes live (non un raster congelato)."""
        if not self.spectra:
            messagebox.showwarning("Save Figure", "Nessuno spettro caricato.")
            return
        path = filedialog.asksaveasfilename(
            initialdir=self.current_dir, defaultextension='.fig.pickle',
            filetypes=[("Matplotlib Figure (pickle)", "*.pickle *.pkl"), ("All Files", "*.*")])
        if not path:
            return
        try:
            self._pulisci_cursore()   # rimuove l'overlay del cursore prima di copiare
            # Salva una copia dimensionata allo stile Origin (single 4:3), non la
            # figura live stirata sul pannello.
            fig_out = pickle.loads(pickle.dumps(self.fig))
            if origin_style is not None and fig_out.axes:
                origin_style.applica_stile_origin(fig_out.axes[0], fig_out, set_size=True, preset='single')
            with open(path, 'wb') as f:
                pickle.dump(fig_out, f)
        except Exception as e:
            traceback.print_exc()
            messagebox.showerror("Save Figure", f"Salvataggio fallito:\n{e}")
            return
        messagebox.showinfo("Save Figure", f"Figura salvata:\n{path}")

    def _carica_plot_editor(self):
        """Importa (una sola volta) il modulo plot_editor: prima un'eventuale copia
        locale in questa cartella, poi il repo fratello PlotStyleKit."""
        if self._pe_module is None:
            import importlib.util
            here = os.path.dirname(os.path.abspath(__file__))
            candidates = [os.path.join(here, 'plot_editor.pyw'),
                          os.path.join(here, '..', 'PlotStyleKit', 'plot_editor.pyw')]
            path = next((p for p in candidates if os.path.isfile(p)), None)
            if path is None:
                raise FileNotFoundError(
                    "plot_editor.pyw non trovato (repo PlotStyleKit mancante accanto a questo progetto)")
            spec = importlib.util.spec_from_file_location('plot_editor', path)
            mod = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(mod)
            self._pe_module = mod
        return self._pe_module

    def apri_editor_figura(self):
        """Apre il Plot Editor direttamente sulla figura corrente, in una nuova
        finestra. L'editor lavora su una copia indipendente (round-trip pickle in
        memoria): così può salvarla/esportarla senza interferire con la vista live,
        che verrebbe comunque ricostruita a ogni aggiornamento."""
        if not self.spectra:
            messagebox.showwarning("Edit Figure", "Nessuno spettro caricato.")
            return
        try:
            pe = self._carica_plot_editor()
        except Exception as e:
            traceback.print_exc()
            messagebox.showerror("Edit Figure", f"plot_editor.pyw non disponibile:\n{e}")
            return
        try:
            self._pulisci_cursore()   # rimuove l'overlay del cursore prima di copiare
            fig_copy = pickle.loads(pickle.dumps(self.fig))
        except Exception as e:
            traceback.print_exc()
            messagebox.showerror("Edit Figure", f"Impossibile duplicare la figura:\n{e}")
            return
        # La figura live è stirata sul pannello: consegna all'editor una copia già
        # dimensionata allo stile Origin (single 4:3) invece che alla dimensione del pannello.
        if origin_style is not None and fig_copy.axes:
            origin_style.applica_stile_origin(fig_copy.axes[0], fig_copy, set_size=True, preset='single')
        top = tk.Toplevel(self.root)
        top.geometry("1300x820")
        editor = pe.PlotEditor(top)
        editor.carica_figura(fig_copy, title="figura corrente")

if __name__ == "__main__":
    root = TkinterDnD.Tk() if HAS_DND else tk.Tk()
    app = LabSpectrumManager(root)
    root.mainloop()
