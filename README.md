# LabSpectrumManager

Viewer and editor for UV-Vis, FTIR and fluorescence spectra: load, overlay and compare multiple
spectra together, with an interactive cursor on the plot and a set of processing tools —
scattering correction (Rayleigh/Mie), adaptive baseline, scaled subtraction between spectra,
averaging, smoothing (boxcar or Savitzky-Golay) and multi-peak deconvolution
(Lorentzian/Gaussian/pseudo-Voigt).

![LabSpectrumManager screenshot](screenshot.png)

## Supported formats

| Extension | Type | Notes |
|---|---|---|
| `.dsp` | UV-Vis (binary) | native parser |
| `.sp` | FTIR PerkinElmer (binary) | native parser (no external tool required) |
| `.csv` | UV-Vis or FTIR | type auto-detected from the X-axis values |

Support for other instrument/vendor-specific formats can be added on request — open an issue with
a sample file.

## Dependencies

`numpy`, `pandas`, `matplotlib`, `tkinter` (standard library); `tkinterdnd2` is optional and
enables drag-and-drop file loading; `scipy` is optional and enables the constrained refinement
step in the scattering correction auto-fit (a plain least-squares fallback is used without it).
[PlotStyleKit](https://github.com/milanone/PlotStyleKit) is an optional sibling repo (cloned as
`../PlotStyleKit` next to this project) providing the shared Origin-like plot style and the
standalone figure editor behind `File → Save Figure (pickle)` / `File → Edit Figure...`; without
it those two menu items are unavailable and the app uses default matplotlib styling.

## Running

```
pythonw LabSpectrumManager.pyw [spectrum_file]
```

## Interface

Three panels side by side: a data table of the merged spectra on the left, a Matplotlib graph
with an interactive crosshair cursor in the center, and a list of loaded spectra (multi-select)
with metadata on the right — which scrolls vertically, since some tool panels (Deconvolution)
can need more room than the window has. The X/Y readout that follows the cursor picks its corner
automatically to avoid the plotted curves, and can be dragged anywhere on the graph.

The data table's columns (one per spectrum, plus X) can be resized by dragging their border, and
the table scrolls both ways once there are more columns or rows than fit. Select rows (click,
Shift/Ctrl to extend, Ctrl+A for all) and Ctrl+C — or right-click — to copy them, header included,
tab-separated, ready to paste into Excel/Origin.

Any spectrum can be hidden from the graph (Hide Selected / Show Selected, or right-click in the
list) without removing it — it stays in the list (greyed out), data table and metadata, just off
the plot. Closing the app (window close button or `File → Exit`) asks for confirmation if there
are calculated results — averages, fits, corrections, trims, pasted data — that haven't been
exported yet with `File → Export CSV`.

`File → Save Figure (pickle)` saves the current graph as a live, re-editable `Figure` object
(not a raster image); `File → Edit Figure...` opens it directly in PlotStyleKit's figure editor
for titles, axis labels, per-line styling, legend placement and publication-size PNG/SVG/PDF
export — see Dependencies.

## Processing tools

- **Scattering correction** — estimates and subtracts a baseline of the form
  `VS · (λ/1000)^-P − M·λ + Offset`, with a constrained auto-fit (SLSQP) over the whole spectrum
  by default: the power term models Rayleigh (P≈4, small particles) to Mie (P down to 1, large
  particles) scattering, the linear term covers residual baseline drift unrelated to scattering
  and is fixed at zero unless unlocked. Theory and operating procedure in
  [`scattering_correction_theory.md`](scattering_correction_theory.md).
- **Linear baseline** (FTIR) — a two-anchor baseline, draggable on the graph, for FTIR spectra.
- **Adaptive baseline** — a lower envelope via iterated moving averages with a minimum
  constraint, adjustable window via a slider (from "flat" to "hugging the signal").
- **Scaled subtraction** — subtracts a reference spectrum B from a spectrum A with an adjustable
  coefficient `k` (`Result = A − k × B`), useful for removing a solvent/blank contribution.
- **Average** of several selected spectra.
- **Smoothing** — repeated moving average (boxcar, ~Gaussian) or Savitzky-Golay filter (local
  polynomial fit, better preserves peak height and shape).
- **Derivative** — 1st or 2nd derivative via Savitzky-Golay (exact derivative of a local
  polynomial fit, in real X units), with the same adjustable smoothing/polynomial-order controls
  as Smoothing, useful for resolving overlapping bands or removing a sloping baseline.
- **Deconvolution** — interactive multi-peak fit (click to add/remove a peak) with a
  Lorentzian, Gaussian or pseudo-Voigt profile; the fit is restricted to a window settable by
  dragging its two boundary lines or typing exact values, same as Trim below.
- **Trim** — crops a spectrum to a window (drag the two vertical lines, or type exact start/end
  values) and stores the selected portion as a new, independent spectrum, leaving the source
  trace untouched; the window edges always snap to a real data point.

## Structure

- `LabSpectrumManager.pyw` — main application (a single `LabSpectrumManager` class)
- `old version and side projects/` — earlier versions (v0-v2) and side projects
  (`scattering.pyw`, `spc_plotter.pyw`) kept as historical reference
- `test_sp_reader.py` — test for the native `.sp` parser (FTIR PerkinElmer)
