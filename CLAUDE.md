# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Running the Application

```bash
# Headless / no console window
pythonw LabSpectrumManager.pyw [optional_spectrum_file]
```

## Dependencies

No requirements file exists. The application relies on these packages being installed in the system Python environment:
- `numpy`, `pandas`, `matplotlib` — data processing and plotting
- `tkinter` — standard library GUI
- `tkinterdnd2` — optional, enables drag-and-drop file loading
- `scipy` — optional, enables the constrained-refinement step in the scattering correction
  auto-fit (`scipy.optimize.minimize`, SLSQP); without it, the auto-fit falls back to the plain
  least-squares solution with a per-point offset clamp
- [PlotStyleKit](https://github.com/milanone/PlotStyleKit) — optional sibling repo (not a pip package). Provides the shared
  "Origin-like" matplotlib style (`origin_style.py`), the standalone figure editor (`plot_editor.pyw`, class `PlotEditor`)
  behind `File → Save Figure (pickle)` / `File → Edit Figure...`, and `plotstylekit.py`, the integration module this
  program talks to. It is found by the common bootstrap block near the top of `LabSpectrumManager.pyw` (the same text as in
  PlotStyleKit's README: change it there, not here): `psk` is that module or None, `PSK_PROBLEM` says why. If it is not found and
  the user has not declined, `_offer_plotstylekit()` asks for the folder 400 ms after the window opens, with a Cancel that is
  remembered; File > PlotStyleKit... (`plotstylekit_dialog()`) changes it later and redraws. Fallbacks without it: default
  matplotlib styling, `salva_figura_pickle()` writes a plain pickle of the figure, and `apri_editor_figura()` offers to
  locate PlotStyleKit (`_need_plotstylekit`). `tests/test_plotstylekit.py` drives all of this with answers in place of the
  dialogs, on a synthetic CSV spectrum. (The block and the methods around it are English, as in every program that uses
  the module, even though this project's own identifiers and comments are Italian.)

No external binaries required — `.sp` (FTIR PerkinElmer) files are parsed natively (see `leggi_sp()` below).

## Architecture

The entire application is a single file (`LabSpectrumManager.pyw`) with one class: `LabSpectrumManager`.

### Data Model

Loaded spectra are stored in `self.spectra` (dict):
```python
{
  'spectrum_name': {
    'df': pd.DataFrame,   # indexed by wavelength/wavenumber, one column per sample
    'info': {'Sample': str, 'Date': str, 'Type': str}
  }
}
```

### UI Layout (3-pane `PanedWindow`)

| Panel | Contents |
|-------|----------|
| Left (400px) | Data table — `ttk.Treeview` grid of merged spectra (X + one column per sample), user-resizable columns, vertical + horizontal scrollbars |
| Center (700px) | Matplotlib graph with interactive cursor overlay |
| Right (400px) | Spectrum listbox (multi-select), remove/clear/hide/show buttons, metadata display, all inside a scrollable canvas (`self.f_right_inner`) since some operation panels (Deconvolution) can need more vertical space than the window has |

### File Format Support

| Extension | Type | Method |
|-----------|------|--------|
| `.dsp` | UV-Vis binary | `leggi_dsp()` — parses binary: wavelength range + point count + Y values |
| `.dsw` | UV-Vis (Varian Cary WinUV) | `leggi_dsw()` — MFC-serialized `TContinuumStore`; data = consecutive float32 (X, Y) pairs located by searching the Start wavelength and validated against the Data Interval; point count from the 'Parameter List' text. Verified on one Cary 500 absorbance scan (Y matches the stored peak labels); first trace only |
| `.sp` | FTIR PerkinElmer | `leggi_sp()` — native binary parser (magic `PEPE`, `DSet2DC1DI` blocks), no external tool |
| `.spc` | Fluorescence (Shimadzu RF-5301) or UV-Vis (Shimadzu UV-1900/UVProbe) | `leggi_spc()` — flat binary for RF-5301; if the file starts with the OLE2 signature `D0CF11E0`, dispatches to `leggi_spc_uvprobe()` (uses the native mini OLE reader `_ole_streams()`; X/Y are float64 streams, X stored descending and re-sorted). UV-1900 verified on one absorbance file only |
| `.csv` | UV-Vis or FTIR | `leggi_csv()` — auto-detects type by threshold: max value > 2000 → FTIR wavenumbers |
| `.csv` | UV-Vis (BMG plate reader) | `leggi_csv()` dispatches to `_leggi_csv_bmg()` when `Test run no.` is in the first 4096 chars — one spectrum per well, optional mean-blank subtraction (asks via `askyesno` unless `sottrai_bianco` is passed); format and limits documented in README "BMG plate-reader CSV"; tests in `tests/test_bmg_reader.py` |

### Key Methods

- `processa_file()` — dispatches to the correct parser based on file extension
- `_popola_tabella_dati()` — rebuilds the Data Table `ttk.Treeview` columns (X + one per
  spectrum) and rows from the merged, rounded DataFrame; called by `aggiorna_vista()`. NaN
  (non-overlapping X ranges between spectra of different resolution/span) renders as an empty
  cell, matching the old CSV-text `to_csv(na_rep='')` behavior it replaced
- `_copia_tabella_dati()` / `_seleziona_tutto_tabella()` — Ctrl+C / Ctrl+A on the Data Table (also
  in its right-click menu) copy the selected rows (header included, tab-separated) to the
  clipboard, replacing the old Text widget's native select-all-and-copy now that it's a Treeview.
  Selecting nothing and pressing Ctrl+C is a no-op (leaves the clipboard untouched)
- Scattering equation: `SCATTER_FORM` (the model as text) and `scatter_equation(VS, P, M, OFF)` (the same with the parameters,
  signs folded into the operators) are module-level and the single place where the formula is written as text; the panel
  shows them (`_sc_eq_label`, refreshed by `_sc_aggiorna_preview()`), `_sc_auto()` adds the window, point count and rms residual
  (`_sc_fit_label`), `_sc_copia_equazione()` fills the clipboard and `_sc_applica()` stores `info['Scatter equation']` next to
  the older low-precision `Scatter VS/P/Slope/Offset` keys. If the model in `_sc_aggiorna_preview()` changes, change the text
  too: `tests/test_scatter_equation.py` evaluates the text against the model
- `LabSpectrumManager.__init__` binds `<ButtonPress-1>` on the `Button` class to `focus_set()`: a Tk button does not take the focus
  when clicked, so a number typed in an entry and not confirmed with Enter (Trim limits, fit range, baseline anchors, ...) was
  not committed before the button's command ran. `tests/test_entry_commit_on_button.py` clicks Apply with real Tk events
- `apri_derivata()` / `_dv_derivata()` — 1st/2nd derivative via `_savgol(..., deriv=1|2, dx=...)`
  (see below); `dx` is `median(abs(diff(x)))`, always positive regardless of whether the spectrum's
  index happens to run ascending or descending, so the derivative's sign follows increasing X by
  convention. Same window/poly controls as Smoothing, but the "Smooth" slider defaults to 0.6 (not
  0.0): a derivative taken with the minimal window is almost pure noise, unlike a lightly-smoothed
  spectrum, so starting near "no smoothing" like Smoothing does would make the first thing the user
  sees look broken
- `apri_shift_scale()` / `_ss_applica()` — Shift / Scale panel: `Y_new = Y * Factor + Offset` on all
  selected spectra (same type required), live preview, result saved as new `<name>_ss` spectra with
  `Factor`/`Offset` in the metadata. The multiplication happens before the offset is added, so the
  offset is in the final Y units
- `_listbox_tasto_destro()` also offers Rename (`_rinomina_spettro()`, one spectrum) and Move up/down/top/
  bottom (`_sposta_spettri()`). The order of `self.spectra` is the order of the list, graph colours and
  CSV export columns, so reordering rebuilds the dict; renaming rebuilds it too (keeps position) and
  must also rename the DataFrame column, `info['Sample']` (if equal to the old name) and the
  `_hidden_spectra` entry. Both refuse to run while an operation panel is open, because panels
  hold spectrum names
- `aggiorna_vista()` — refreshes all three panels after any data change; auto-selects the
  listbox entry when exactly one spectrum is loaded, so single-selection tools (Trim,
  Smoothing, Deconvolution, ...) work without an explicit click first
- `_evidenzia_selezionati()` — on `<<ListboxSelect>>` (and at the end of `aggiorna_vista()`, since `selection_set()` fires no event) thickens the plot lines of the selected spectra (lw 2.6 vs 1.5) without redrawing everything. Only lines tagged `_lsm_main` by `aggiorna_vista()` are touched, so traces drawn by operation panels keep their own widths; legend handles keep the default width until the next redraw
- `on_mouse_move()` / `on_mouse_leave()` — interactive crosshair cursor tracking on plot
- `carica_da_dialog()` / `handle_drop()` — file import entry points (dialog and DnD)
- `esporta_csv()` — exports merged spectra to CSV
- `salva_figura_pickle()` / `apri_editor_figura()` — save the current figure as a live
  `Figure` pickle / open it in PlotStyleKit's `PlotEditor` (see Dependencies above) through `psk.save_pickle_dialog()` /
  `psk.edit_dialog()`; both work on an in-memory pickle round-trip copy, restyled to the Origin `single` preset (`unpin=False`:
  this view does not pin its legend), so the live panel view is never mutated
- `apri_trim()` / `_tr_applica()` — crop a spectrum to a window (draggable vertical lines or
  typed values in the entry boxes, both snapped to the nearest real data point — never an
  interpolated x) and store the selection as a new, independent spectrum; the source trace is
  left untouched
- `hide_selected()` / `show_selected()` (and the listbox's right-click Hide/Show) — exclude a
  spectrum from the graph and cursor readout without removing it; still listed (greyed out) and
  still present in the data table and metadata panel
- `apri_deconvoluzione()` — the fit-restriction window (`self._dc_range`, two draggable vertical
  lines) has entry boxes like Trim's, same real-data-point snapping. Left-click both adds a peak
  and starts a line drag if near a line; `_dc_drag_stop()` only counts it as a drag (and skips
  adding a peak) if the mouse actually moved (`self._dc_drag_moved`) — otherwise a stationary
  click near either line is treated as a normal add-peak click. Without this, a peak near either
  edge of the range is unclickable, which matters now that the range defaults to the whole
  spectrum (both lines sit right at the plot's outer edges)
- `_dc_ridisegna()` — the drawn fit/guess curve and each peak component are masked to
  `self._dc_range` (`curve_fit` in `_dc_fit()` already only used that window); only the raw
  spectrum trace is drawn across the full range, for context
- `_dc_add_mode` (checkbox "Add/remove peaks") gates click-to-add and right-click-remove
  in `_dc_press()`/`_dc_drag_stop()` — the range-line drag stays active either way. Toggling it
  on, and opening the panel at all, call `_dc_disattiva_toolbar_mode()`: matplotlib's pan/zoom
  toolbar mode silently swallows canvas clicks while active, which is a common cause of "clicking
  does nothing" after zooming in to inspect a peak before adding it
- The plot canvas's `<Button-3>` binding (`plot_widget.bind('<Button-3>', self._menu_grafico,
  add='+')`) MUST keep `add='+'`. `FigureCanvasTk.__init__` already binds `<Button-1/2/3>` on the
  same Tk widget to generate matplotlib's `button_press_event`; binding `<Button-3>` again
  without `add='+'` replaces that binding instead of adding to it, so no single right-click ever
  reaches any `mpl_connect('button_press_event', ...)` handler again — e.g. deconvolution's
  right-click-to-remove — while `<Double-Button-3>` (bound separately by matplotlib) still works,
  which is why the symptom looks like "only double right-click removes a peak"
- `on_exit()` — wired to both `File → Exit` and the window close button; asks for confirmation
  when `self._dirty` is set. `self._dirty` is set by every operation that creates a derived or
  clipboard-pasted spectrum (average, scattering correction, subtraction, normalize, baseline,
  smoothing, deconvolution, trim, paste), and cleared by `esporta_csv()` and `clear_all()`

Method names follow Italian conventions (`leggi` = read, `aggiorna` = update, `carica` = load).

## Editing conventions
- Edits must be surgical and non-destructive
- Never refactor or rename existing methods unless explicitly asked
- Never change the public interface without explicit instruction
- Preserve all existing comments and docstrings
- When in doubt, ask before modifying
- User-facing widget text — button/checkbox/label captions, panel titles — is English (the exit
  confirmation, panel titles like "Scattering Correction"/"Trim", "Fit range:", "From"/"To", etc.
  all are). `messagebox` body text stays Italian (the app's established mixed convention: English
  chrome, Italian prose). Code comments and identifiers stay Italian either way.
- Any new widget that belongs in the right panel (a new operation button, a new panel frame) must
  be parented to `self.f_right_inner`, not `self.f_right` — the latter is now just the outer
  frame holding the scrollable canvas + scrollbar; parenting to it directly would place the
  widget outside the scrollable area.
