# Teoria della correzione scattering

La correzione rimuove dallo spettro misurato un **baseline artificiale** generata da fenomeni di scattering della luce, che non appartiene alla vera risposta ottica del campione.

---

## Origine fisica

Quando la luce attraversa una soluzione contenente particelle (aggregati, colloidi, polveri in sospensione, membrane), una frazione viene **deviata** invece di essere assorbita. Lo spettrofotometro registra questa deviazione come "assorbanza apparente", sovrapposta al segnale reale.

---

## Legge di potenza (Rayleigh / Mie)

L'intensità dello scattering dipende dalla lunghezza d'onda secondo una legge di potenza:

$$S(\lambda) = A_{eff} \cdot \lambda^{-P}$$

| P | Regime | Particelle tipiche |
|---|--------|-------------------|
| **4** | Rayleigh | Molto più piccole della λ (proteine, nanoparticelle < 50 nm). Decade rapidamente verso il rosso. |
| **2–3** | Mie | Comparabili alla λ (liposomi, membrane, aggregati). Più piatto e persistente. |
| **< 2** | Mie estremo / scattering multiplo | Particelle grandi. |

In pratica P si sceglie empiricamente osservando quanto è "ripida" la coda verso il blu.

---

## Componente lineare

In molti campioni reali si sovrappone una **deriva lineare di baseline** indipendente dallo scattering classico, dovuta a:

- disuniformità ottiche della cuvetta
- differenze di indice di rifrazione tra campione e riferimento
- scattering multiplo a bassa frequenza

Questa componente non segue una legge di potenza e non viene rimossa dal termine $A_{eff} \cdot \lambda^{-P}$. Va trattata separatamente come termine lineare $m \cdot \lambda$.

---

## Formula completa implementata

$$\text{correction}(\lambda) = A \cdot 10^{VS} \cdot \left(\frac{\lambda}{1000}\right)^{-P} + m \cdot \lambda + \text{Offset}$$

| Parametro | Ruolo | Quando usarlo |
|---|---|---|
| **A** | Ampiezza dello scattering | Sempre, se c'è scattering |
| **VS** | Fattore di scala logaritmico per A | Quando A è molto piccolo o molto grande |
| **P** | Esponente (2, 3, 4) | Dipende dalla dimensione delle particelle |
| **m** | Pendenza della componente lineare | Se la baseline è inclinata dopo il termine di potenza |
| **Offset** | Traslazione verticale costante | Per azzerare la baseline a λ alte dove lo scattering è trascurabile |

---

## Procedura operativa consigliata

1. Partire da **P = 4** (Rayleigh) e aumentare A fino a correggere la coda UV
2. Se residua una pendenza, abbassare P a 3 o 2 (Mie)
3. Se rimane un tilt, aggiustare **m** (Slope)
4. Usare **Offset** per traslare la baseline a zero in una zona priva di assorbanza reale (tipicamente > 750 nm per la maggior parte dei cromofori biologici)

---

> La stima dei parametri è visiva/empirica: non esiste un metodo automatico universale perché la proporzione tra scattering Rayleigh, Mie e deriva lineare dipende dal campione specifico.
