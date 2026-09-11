# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Running the Application

```bash
# Via Windows batch launcher (passes optional file path arg)
LabSpectrumManager.bat [optional_spectrum_file]

# Direct Python (headless/no console window)
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
- [PlotStyleKit](https://github.com/milanone/PlotStyleKit) — optional sibling repo (not a pip
  package), expected as `../PlotStyleKit` next to this project's folder. Provides the shared
  "Origin-like" matplotlib style (`origin_style.py`) and the standalone figure editor
  (`plot_editor.pyw`, class `PlotEditor`) behind `File → Save Figure (pickle)` / `File → Edit
  Figure...`. Loaded by path at runtime (local copy in this folder first, then the sibling
  folder); without it the app runs normally with default matplotlib styling and those two menu
  items unavailable — see the one-time startup warning in `__init__`. Same loading pattern and
  reference implementation as `KleistekManager.pyw`.

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
| `.sp` | FTIR PerkinElmer | `leggi_sp()` — native binary parser (magic `PEPE`, `DSet2DC1DI` blocks), no external tool |
| `.csv` | UV-Vis or FTIR | `leggi_csv()` — auto-detects type by threshold: max value > 2000 → FTIR wavenumbers |

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
- `apri_derivata()` / `_dv_derivata()` — 1st/2nd derivative via `_savgol(..., deriv=1|2, dx=...)`
  (see below); `dx` is `median(abs(diff(x)))`, always positive regardless of whether the spectrum's
  index happens to run ascending or descending, so the derivative's sign follows increasing X by
  convention. Same window/poly controls as Smoothing, but the "Smooth" slider defaults to 0.6 (not
  0.0): a derivative taken with the minimal window is almost pure noise, unlike a lightly-smoothed
  spectrum, so starting near "no smoothing" like Smoothing does would make the first thing the user
  sees look broken
- `aggiorna_vista()` — refreshes all three panels after any data change; auto-selects the
  listbox entry when exactly one spectrum is loaded, so single-selection tools (Trim,
  Smoothing, Deconvolution, ...) work without an explicit click first
- `on_mouse_move()` / `on_mouse_leave()` — interactive crosshair cursor tracking on plot
- `carica_da_dialog()` / `handle_drop()` — file import entry points (dialog and DnD)
- `esporta_csv()` — exports merged spectra to CSV
- `salva_figura_pickle()` / `apri_editor_figura()` — save the current figure as a live
  `Figure` pickle / open it in PlotStyleKit's `PlotEditor` (see Dependencies above); both work
  on an in-memory pickle round-trip copy, restyled to the Origin `single` preset, so the live
  panel view is never mutated
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
