# Teoria della correzione scattering

La correzione rimuove dallo spettro misurato un **baseline artificiale** generata da fenomeni di scattering della luce, che non appartiene alla vera risposta ottica del campione.

---

## Parte 1 — Teoria fisica

### Origine fisica

Quando la luce attraversa una soluzione contenente particelle (aggregati, colloidi, polveri in sospensione, membrane, cellule), una frazione viene **deviata** invece di essere assorbita. Lo spettrofotometro registra questa deviazione come "assorbanza apparente", sovrapposta al segnale reale.

### Legge di potenza (Rayleigh / Mie)

L'intensità dello scattering dipende dalla lunghezza d'onda secondo una legge di potenza:

$$S(\lambda) \propto \lambda^{-P}$$

| P | Regime | Particelle tipiche |
|---|--------|-------------------|
| **4** | Rayleigh | Molto più piccole della λ (proteine, nanoparticelle < 50 nm). Decade rapidamente verso il rosso. |
| **2–3** | Mie | Comparabili alla λ (liposomi, membrane, aggregati). Più piatto e persistente. |
| **1** | Mie estremo / scattering multiplo | Particelle grandi rispetto a λ (es. cellule intere, ~1–2 µm). Curva quasi piatta. |

In generale P non è un parametro libero su tutto l'asse reale: sotto P≈1 la legge di potenza `λ^-P` si appiattisce verso una costante e smette di essere distinguibile, in un fit, dal semplice offset verticale — il modello perde di significato fisico ancora prima che matematico. Per questo motivo P ha senso solo per P ≥ 1 (si veda la Parte 2 per il vincolo effettivo usato dal programma).

### Deriva residua non riconducibile allo scattering

In alcuni campioni reali può sovrapporsi una **deriva di baseline** che non segue la legge di potenza dello scattering, dovuta ad esempio a:

- disuniformità ottiche della cuvetta
- differenze di indice di rifrazione tra campione e riferimento
- imperfezioni strumentali a bassa frequenza

Quando presente, in prima approssimazione (su un intervallo spettrale non troppo esteso) può essere descritta come una componente lineare in λ, da sottrarre in aggiunta al termine di scattering.

---

## Parte 2 — Implementazione nel programma

### Formula usata

Il pannello "Scattering Correction" calcola e sottrae dallo spettro la seguente correzione:

$$\text{correction}(\lambda) = VS \cdot \left(\frac{\lambda}{1000}\right)^{-P} - M \cdot \lambda + \text{OFF}$$

| Parametro | Ruolo | Range nel programma |
|---|---|---|
| **VS** | Ampiezza dello scattering (vincolata ≥ 0) | Slider continuo, calibrato sull'ampiezza dello spettro caricato |
| **P** | Esponente della legge di potenza | Slider continuo **1.0 – 5.0**, con preset rapidi 1/2/3/4 |
| **M** | Pendenza della deriva lineare residua | Slider continuo; **0 di default**, bloccato da "fix slope" (vedi sotto) |
| **OFF** | Traslazione verticale costante | Slider continuo, per azzerare la baseline dove lo scattering è trascurabile |

Il limite inferiore P ≥ 1 non è arbitrario: è il punto in cui il termine di potenza degenera (vedi Parte 1). Un test su uno spettro reale (`spettro_sphaer.csv`) ha mostrato che abbassare il limite a 0.1 riduce l'SSE nella finestra di fit solo di ~4‰, ma fa divergere VS e OFF verso valori enormi e di segno opposto che si cancellano a vicenda — un guadagno numerico irrilevante pagato con un fit non identificabile.

### Il termine lineare M: perché esiste ma è bloccato di default

La componente lineare `M·λ` intercetta derive di baseline non spiegabili dal termine di potenza (vedi Parte 1). Con P libero nell'intervallo [1, 5], il termine di potenza da solo copre la quasi totalità dei casi osservati: quando P è stimato correttamente, M resta vicino a zero anche quando è incluso nel fit.

Per questo **M è fissato a 0 di default** (checkbox "fix slope" attivo) e va sbloccato solo se l'auto-fit non converge in modo soddisfacente lasciando P libero — cioè come rete di sicurezza, non come parametro di uso corrente.

### Auto-fit vincolato

Il pulsante "Auto (fit window)" stima i parametri risolvendo un problema di minimi quadrati vincolato (SLSQP, richiede `scipy`) sulla finestra delimitata dalle due linee verticali trascinabili — che per default coprono l'intero spettro: l'algoritmo vincolato converge in modo affidabile su tutta la banda:

- il vincolo di disuguaglianza impone che la traccia di correzione non superi mai lo spettro misurato (traccia ≤ spettro, punto per punto);
- se **"fix P" non è spuntato**, l'esponente P viene ottimizzato insieme a VS, OFF (e M, se sbloccato) entro i limiti [1, 5]; se è spuntato, P resta al valore corrente e il fit si riduce a un problema lineare (risolto con jacobiano analitico, più veloce e stabile);
- se **"fix slope" non è spuntato**, M viene incluso tra i parametri liberi del fit; se è spuntato (default), resta bloccato al valore corrente (0 di norma) e il fit ottimizza solo VS, OFF (ed eventualmente P);
- senza `scipy` installato, l'auto-fit si riduce alla soluzione ai minimi quadrati non vincolata, con un abbassamento dell'offset a posteriori per garantire comunque traccia ≤ spettro.

### Procedura operativa consigliata

1. Selezionare lo spettro e aprire "Scattering Correction": le linee di finestra coprono già tutto lo spettro.
2. Premere **Auto (fit window)** lasciando P libero e "fix slope" attivo (default): nella maggior parte dei casi è sufficiente.
3. Se il residuo mostra ancora una pendenza sistematica, provare a sbloccare **"fix slope"** e rilanciare Auto.
4. Se si vuole forzare un regime fisico noto a priori (es. P=4 per particelle molto piccole), spuntare **"fix P"**, impostare il preset desiderato e rilanciare Auto: verranno ottimizzati solo VS (e OFF, M).
5. Usare **Offset** per rifiniture manuali, tipicamente per azzerare la baseline in una zona priva di assorbanza reale (es. > 750 nm per molti cromofori biologici).

### Offset bloccato a 0 e curva P(d) dei liposomi

Per particelle che non assorbono (liposomi) la legge pura `VS·(λ/1000)^(−P)` è il modello giusto, e un offset libero ne assorbe parte della curvatura abbassando P di 0.3-0.7 (a 60-120 nm). Spuntando **"fix offset = 0"**, **Auto** esegue un fit ai minimi quadrati semplice (senza il vincolo traccia ≤ spettro, che senza offset da abbassare spingerebbe la traccia sotto i dati e distorcerebbe P), con P libero e pendenza bloccata a 0. Il residuo mostra se la baseline è davvero nulla.

Il pulsante **"P vs size (theory)..."** calcola, con la teoria di Mie per una vescicola (guscio di 4.5 nm, Δn 0.12, acqua dentro e fuori) e popolazione log-normale (σ 0.2 e 0.3), il P che lo stesso fit dà sullo spettro teorico di vescicole di ogni diametro, nella finestra di fit corrente, con l'angolo di semiapertura di raccolta dello strumento (10° per lo spettrofotometro usato finora; la luce diffusa entro quell'angolo arriva al rivelatore e non conta come estinzione). La retta rossa è il P del pannello: dove incontra la curva sono le dimensioni compatibili.

- Rayleigh (P = 4) vale solo sotto circa 30 nm; a 120 nm P è ≈ 2.6-2.7, quindi P = 4 non è il valore atteso.
- Tra circa 40 e 130 nm la curva scende in modo netto (≈ −0.5 ogni 20 nm): P determina la dimensione. Sopra 130-150 nm è piatta (2.4-2.6): P non la determina, serve l'ampiezza a concentrazione nota o il DLS.
- La curva dipende da angolo di raccolta, spessore e indice del guscio, ma poco: l'angolo conta soprattutto sopra 200 nm.
- Richiede `scipy` (solo per questo calcolo). Verificato su tre campioni POPC-colesterolo-PEG, confrontati con il DLS: un campione di 177 nm concorda, uno di 126 nm (DLS) dà ~80 nm dallo spettro. Il modello è unilamellare: vescicole multilamellari hanno ampiezza maggiore.

Il pulsante **"Spectrum simulator..."** sovrappone allo spettro misurato lo spettro teorico in scala assoluta, senza ampiezza libera: lo calcola dalla concentrazione di lipide (stock, diluizione, perdite, cammino ottico), dal diametro e dalla larghezza della popolazione, dal numero di lamelle, dal guscio e dall'angolo di raccolta. Sotto compare il residuo (misurato − teoria) con rms e il fattore di ampiezza k che li farebbe coincidere (k = 1: la scala assoluta concorda). **Fit size** cerca il diametro con rms minimo a k dato; **Fit size + amplitude** risolve anche k. La scala assoluta dipende da spessore del doppio strato × Δn (fattore 2-3): conviene tarare k su campioni estrusi unilamellari e usarlo poi per contare le lamelle degli altri. Il modello multilamellare usa lamelle equidistanti (7 nm di ripetizione) e Δn mediato sul guscio.

La composizione (mol% di colesterolo e di lipide-PEG, il resto è POPC) entra nel simulatore così: massa molare media (es. 825 g/mol per 85/10/5 con DSPE-PEG2000), area per molecola ridotta del 9% ogni 10 mol% di colesterolo, spessore del bilayer +0.4 nm ogni 10 mol%, e la corona di PEG come film di polimero puro (volume della catena = massa/1.2 g/cm³, una catena ogni 20 lipidi al 5 mol%, circa 0.23 nm per foglietto) con indice 0.127 sopra l'acqua, che si somma alla polarizzabilità del guscio (la corona è molto più sottile della lunghezza d'onda: conta solo il totale). Sono stime di ordine di grandezza. Con 85/10/5 il fattore k del campione 200 nm scende da 1.8 a 1.5; le dimensioni non cambiano.

> La stima dei parametri combina un auto-fit vincolato, che copre il caso generale, con l'intervento manuale — lock di P e/o M, regolazione fine dei singoli slider — per casi limite o per imporre un regime fisico noto a priori.
