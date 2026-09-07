# LabSpectrumManager

Viewer ed editor per spettri UV-Vis, FTIR e di fluorescenza: carica, sovrappone e
confronta più spettri insieme, con cursore interattivo sul grafico e un set di strumenti di
elaborazione — correzione scattering (Rayleigh/Mie), baseline adattiva, sottrazione scalata tra
spettri, media, smoothing (boxcar o Savitzky-Golay) e deconvoluzione multi-picco
(Lorentz/Gauss/pseudo-Voigt).

## Formati supportati

| Estensione | Tipo | Note |
|---|---|---|
| `.dsp` | UV-Vis (binario) | parser nativo |
| `.sp` | FTIR PerkinElmer (binario) | parser nativo (nessun tool esterno richiesto) |
| `.csv` | UV-Vis o FTIR | il tipo viene rilevato automaticamente dai valori sull'asse X |

Il supporto per altri formati specifici (altri strumenti/produttori) può essere aggiunto su
richiesta — apri una issue con un file di esempio.

## Dipendenze

`numpy`, `pandas`, `matplotlib`, `tkinter` (libreria standard); `tkinterdnd2` opzionale, abilita il
trascinamento dei file nella finestra.

## Avvio

```
pythonw LabSpectrumManager.pyw [file_spettro]
```

## Interfaccia

Tre pannelli affiancati: tabella dati (spettri uniti come testo tab-separated) a sinistra, grafico
Matplotlib con cursore a croce interattivo al centro, elenco spettri caricati (multi-selezione) con
metadati a destra.

## Strumenti di elaborazione

- **Correzione scattering** — stima e sottrae una baseline della forma
  `A · 10^VS · (λ/1000)^-P + m·λ + Offset`: il termine di potenza modella lo scattering Rayleigh
  (P=4, particelle piccole) o Mie (P=2–3, particelle comparabili alla lunghezza d'onda), il
  termine lineare copre derive di baseline non legate allo scattering. Dettagli teorici e
  procedura operativa in [`scattering_correction_theory.md`](scattering_correction_theory.md).
- **Baseline lineare** (FTIR) — baseline a due ancore trascinabili sul grafico, per spettri FTIR.
- **Baseline adattiva** — inviluppo inferiore via medie mobili iterate con vincolo di minimo,
  finestra regolabile con uno slider (da "piatta" a "aderente al segnale").
- **Sottrazione scalata** — sottrae uno spettro di riferimento B da uno spettro A con un
  coefficiente `k` regolabile (`Result = A − k × B`), utile per rimuovere il contributo del
  solvente/bianco.
- **Media** di più spettri selezionati.
- **Smoothing** — media mobile ripetuta (boxcar, ~gaussiana) oppure filtro Savitzky-Golay
  (fit polinomiale locale, preserva meglio altezza e forma dei picchi).
- **Deconvoluzione** — fit multi-picco interattivo (clic per aggiungere/rimuovere un picco) con
  profilo Lorentziano, Gaussiano o pseudo-Voigt.

## Struttura

- `LabSpectrumManager.pyw` — applicazione principale (un'unica classe `LabSpectrumManager`)
- `old version and side projects/` — versioni precedenti (v0-v2) e side-project (`scattering.pyw`,
  `spc_plotter.pyw`) tenuti come riferimento storico
- `test_sp_reader.py` — test del parser nativo per il formato `.sp` (FTIR PerkinElmer)
