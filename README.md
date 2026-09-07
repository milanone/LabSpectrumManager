# LabSpectrumManager

Viewer per spettri di laboratorio (UV-Vis, FTIR, fluorescenza): carica, sovrappone e confronta più
spettri insieme, con cursore interattivo sul grafico e correzione dello scattering (Rayleigh/Mie)
per rimuovere la baseline artificiale introdotta da particelle in sospensione (liposomi, membrane,
aggregati) nei campioni UV-Vis.

## Formati supportati

| Estensione | Tipo | Note |
|---|---|---|
| `.dsp` | UV-Vis (binario) | parser nativo |
| `.sp` | FTIR PerkinElmer (binario) | parser nativo (nessun tool esterno richiesto) |
| `.csv` | UV-Vis o FTIR | il tipo viene rilevato automaticamente dai valori sull'asse X |

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

## Correzione scattering

Il pannello dedicato stima e sottrae una baseline della forma
`A · 10^VS · (λ/1000)^-P + m·λ + Offset` — il termine di potenza modella lo scattering
Rayleigh (P=4, particelle piccole) o Mie (P=2–3, particelle comparabili alla lunghezza d'onda),
il termine lineare copre derive di baseline non legate allo scattering. Dettagli teorici e
procedura operativa in [`scattering_correction_theory.md`](scattering_correction_theory.md).

## Struttura

- `LabSpectrumManager.pyw` — applicazione principale (un'unica classe `LabSpectrumManager`)
- `old version and side projects/` — versioni precedenti (v0-v2) e side-project (`scattering.pyw`,
  `spc_plotter.pyw`) tenuti come riferimento storico
- `test_sp_reader.py` — test del parser nativo per il formato `.sp` (FTIR PerkinElmer)
