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

Method names follow Italian conventions (`leggi` = read, `aggiorna` = update, `carica` = load).

## Editing conventions
- Edits must be surgical and non-destructive
- Never refactor or rename existing methods unless explicitly asked
- Never change the public interface without explicit instruction
- Preserve all existing comments and docstrings
- When in doubt, ask before modifying
