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
| Left (400px) | Data table — merged spectra as tab-separated text |
| Center (700px) | Matplotlib graph with interactive cursor overlay |
| Right (400px) | Spectrum listbox (multi-select), remove/clear buttons, metadata display |

### File Format Support

| Extension | Type | Method |
|-----------|------|--------|
| `.dsp` | UV-Vis binary | `leggi_dsp()` — parses binary: wavelength range + point count + Y values |
| `.sp` | FTIR PerkinElmer | `leggi_sp()` — native binary parser (magic `PEPE`, `DSet2DC1DI` blocks), no external tool |
| `.csv` | UV-Vis or FTIR | `leggi_csv()` — auto-detects type by threshold: max value > 2000 → FTIR wavenumbers |

### Key Methods

- `processa_file()` — dispatches to the correct parser based on file extension
- `aggiorna_vista()` — refreshes all three panels after any data change
- `on_mouse_move()` / `on_mouse_leave()` — interactive crosshair cursor tracking on plot
- `carica_da_dialog()` / `handle_drop()` — file import entry points (dialog and DnD)
- `esporta_csv()` — exports merged spectra to CSV
- `salva_figura_pickle()` / `apri_editor_figura()` — save the current figure as a live
  `Figure` pickle / open it in PlotStyleKit's `PlotEditor` (see Dependencies above); both work
  on an in-memory pickle round-trip copy, restyled to the Origin `single` preset, so the live
  panel view is never mutated

Method names follow Italian conventions (`leggi` = read, `aggiorna` = update, `carica` = load).

## Editing conventions
- Edits must be surgical and non-destructive
- Never refactor or rename existing methods unless explicitly asked
- Never change the public interface without explicit instruction
- Preserve all existing comments and docstrings
- When in doubt, ask before modifying
