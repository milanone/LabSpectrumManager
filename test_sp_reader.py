"""
test_sp_reader.py
-----------------
Parser e visualizzatore standalone per file binari .SP di PerkinElmer FTIR.

Formato reverse-engineered da:
  github.com/kutukvpavel/PerkinElmer2CSV

Struttura file:
  Bytes 0-3   : magic "PEPE"
  Bytes 4-43  : descrizione ASCII (40 byte)
  Bytes 44+   : sequenza di blocchi top-level

Ogni blocco top-level:
  Int16  id
  Int32  length
  bytes[length]  data

Il blocco con id=120 (DSet2DC1DI) contiene i dati spettrali
come sequenza di member-block:
  Int16  id
  Int32  length
  Int16  typecode
  bytes[length-2]  data

Member id rilevanti (little-endian int16, valori negativi):
  -29838  DataSetAbscissaRange  → 2 double: StartX, EndX
  -29835  DataSetNumPoints      → 1 int32 : numero punti
  -29828  DataSetData           → int32 count + count*double: valori Y

Uso:
  python test_sp_reader.py file.sp [file2.sp ...]
  python test_sp_reader.py          (apre dialog di selezione file)
"""

import struct
import sys
import os
import re
import datetime
import numpy as np
import matplotlib.pyplot as plt
import tkinter as tk
from tkinter import filedialog


# ---------------------------------------------------------------------------
# Costanti del formato
# ---------------------------------------------------------------------------
MAGIC           = b'PEPE'
HEADER_SIZE     = 44          # 4 (magic) + 40 (descrizione)
DSET2DC1DI      = 120         # id blocco wrapper dati spettrali
ABSCISSA_RANGE  = -29838      # DataSetAbscissaRange
NUM_POINTS      = -29835      # DataSetNumPoints
DATA            = -29828      # DataSetData


# ---------------------------------------------------------------------------
def parse_sp(path):
    """
    Legge un file .SP di PerkinElmer e restituisce (x_vals, y_vals, description).
    x_vals e y_vals sono numpy array; description è una stringa.
    Lancia ValueError se il file non è riconosciuto o i dati mancano.
    """
    with open(path, 'rb') as f:
        raw = f.read()

    if len(raw) < HEADER_SIZE or raw[:4] != MAGIC:
        raise ValueError(f"'{os.path.basename(path)}': magic 'PEPE' mancante — non è un file SP valido")

    description = raw[4:44].decode('ascii', errors='replace').rstrip('\x00').strip()

    x_start = x_end = n_points = None
    y_data  = None
    date_str = None

    pos = HEADER_SIZE
    while pos + 6 <= len(raw):
        block_id  = struct.unpack_from('<h', raw, pos)[0]
        block_len = struct.unpack_from('<i', raw, pos + 2)[0]
        pos += 6
        block_data = raw[pos:pos + block_len]
        pos += block_len

        if block_id != DSET2DC1DI:
            continue  # salta blocchi che non sono dati spettrali

        # --- parsifica member-block interni ---
        mpos = 0
        while mpos + 8 <= len(block_data):
            m_id  = struct.unpack_from('<h', block_data, mpos)[0]
            m_len = struct.unpack_from('<i', block_data, mpos + 2)[0]
            mpos += 6
            # primi 2 byte = TypeCode, il resto = dati effettivi
            m_data = block_data[mpos + 2 : mpos + m_len]
            mpos  += m_len

            if m_id == -29825:  # DataSetHistoryRecord
                date_str = _extract_date(m_data)

            if m_id == ABSCISSA_RANGE and len(m_data) >= 16:
                x_start = struct.unpack_from('<d', m_data, 0)[0]
                x_end   = struct.unpack_from('<d', m_data, 8)[0]
                print(f"  AbscissaRange: {x_start:.4f} → {x_end:.4f} cm⁻¹")

            elif m_id == NUM_POINTS and len(m_data) >= 4:
                n_points = struct.unpack_from('<i', m_data, 0)[0]
                print(f"  NumPoints: {n_points}")

            elif m_id == DATA and len(m_data) >= 4:
                # Il primo Int32 è la dimensione in BYTE dell'array dati,
                # non il numero di elementi. Dividendo per 8 si ottiene il
                # numero di double (little-endian, 64-bit).
                byte_len = struct.unpack_from('<i', m_data, 0)[0]
                n_elem   = byte_len // 8
                if len(m_data) < 4 + n_elem * 8:
                    raise ValueError(f"DataSetData: buffer atteso {4 + n_elem*8}B, disponibile {len(m_data)}B")
                y_data = np.array(struct.unpack_from(f'<{n_elem}d', m_data, 4))
                print(f"  DataSetData: {n_elem} double letti ({byte_len} byte)")

        break  # trovato il blocco dati, non serve proseguire

    if x_start is None or x_end is None:
        raise ValueError("AbscissaRange non trovato nel file")
    if y_data is None:
        raise ValueError("DataSetData non trovato nel file")

    if n_points is None:
        n_points = len(y_data)

    x_vals = np.linspace(x_start, x_end, n_points)
    return x_vals, y_data, description, date_str


def _extract_date(history_data):
    """Estrae la data di misura originale dal DataSetHistoryRecord.

    Cerca la voce 'Created by Instrument' e il timestamp ctime che la segue
    (formato: 'ddd mmm dd HH:MM:SS YYYY'), restituisce 'YYYY-MM-DD HH:MM'.
    """
    try:
        text = history_data.decode('latin-1', errors='replace')
        # Il timestamp segue 'Created by Instrument' entro ~40 caratteri
        m = re.search(
            r'Created by Instrument.{1,10}([A-Za-z]{3} [A-Za-z]{3} {1,2}\d{1,2} \d{2}:\d{2}:\d{2} \d{4})',
            text, re.DOTALL
        )
        if m:
            ts = m.group(1).strip()
            dt = datetime.datetime.strptime(ts, '%a %b %d %H:%M:%S %Y')
            return dt.strftime('%Y-%m-%d %H:%M')
    except Exception:
        pass
    return None


# ---------------------------------------------------------------------------
# Nomi simbolici dei blocchi top-level e dei member (dal sorgente C#)
BLOCK_NAMES = {
    120: 'DSet2DC1DI', 121: 'HistoryRecord', 122: 'InstrHdrHistoryRecord',
    123: 'InstrumentHeader', 124: 'IRInstrumentHeader',
    125: 'UVInstrumentHeader', 126: 'FLInstrumentHeader',
}
MEMBER_NAMES = {
    -29839: 'DataSetDataType',      -29838: 'DataSetAbscissaRange',
    -29837: 'DataSetOrdinateRange', -29836: 'DataSetInterval',
    -29835: 'DataSetNumPoints',     -29834: 'DataSetSamplingMethod',
    -29833: 'DataSetXAxisLabel',    -29832: 'DataSetYAxisLabel',
    -29831: 'DataSetXAxisUnitType', -29830: 'DataSetYAxisUnitType',
    -29829: 'DataSetFileType',      -29828: 'DataSetData',
    -29827: 'DataSetName',          -29826: 'DataSetChecksum',
    -29825: 'DataSetHistoryRecord', -29824: 'DataSetInvalidRegion',
    -29823: 'DataSetAlias',         -29822: 'DataSetVXIRAccyHdr',
    -29821: 'DataSetVXIRQualHdr',   -29820: 'DataSetEventMarkers',
}

def _try_decode(data):
    """Prova a estrarre testo leggibile dai byte del member."""
    # Formato stringa PE: Int16 length + ASCII
    if len(data) >= 2:
        slen = struct.unpack_from('<h', data, 0)[0]
        if 0 < slen <= len(data) - 2:
            candidate = data[2:2 + slen].decode('ascii', errors='replace')
            if candidate.isprintable():
                return repr(candidate)
    # Fallback: decodifica raw
    try:
        s = data.decode('latin-1').rstrip('\x00').strip()
        if s.isprintable() and len(s) > 1:
            return repr(s)
    except Exception:
        pass
    return None

def decode_history(data):
    """Tenta di estrarre data/ora e testo dal member DataSetHistoryRecord."""
    print(f"    --- DataSetHistoryRecord ({len(data)} byte) ---")
    # Dump hex delle prime 32 byte per ispezione
    print(f"    hex[0:32]: {data[:32].hex(' ')}")

    # Tentativo 1: Windows FILETIME (Int64, 100-ns dal 1601-01-01)
    if len(data) >= 8:
        import datetime
        ft = struct.unpack_from('<q', data, 0)[0]
        try:
            dt = datetime.datetime(1601, 1, 1) + datetime.timedelta(microseconds=ft // 10)
            if 1990 < dt.year < 2100:
                print(f"    -> FILETIME[0]: {dt.strftime('%Y-%m-%d %H:%M:%S')}")
        except Exception:
            pass

    # Tentativo 2: FILETIME a offset 4 o 8
    for off in (4, 8, 12):
        if len(data) >= off + 8:
            ft = struct.unpack_from('<q', data, off)[0]
            try:
                import datetime
                dt = datetime.datetime(1601, 1, 1) + datetime.timedelta(microseconds=ft // 10)
                if 1990 < dt.year < 2100:
                    print(f"    -> FILETIME[{off}]: {dt.strftime('%Y-%m-%d %H:%M:%S')}")
            except Exception:
                pass

    # Tentativo 3: cerca testo leggibile (date, nomi) nel blocco
    try:
        text = data.decode('latin-1', errors='replace')
        # Cerca pattern data (es. "2025", "06/07", ...)
        import re
        for m in re.finditer(r'[\x20-\x7e]{4,}', text):
            s = m.group().strip()
            if s:
                print(f"    testo[{m.start()}]: {repr(s)}")
    except Exception:
        pass
    print()


def dump_sp(path):
    """Stampa la struttura completa di un file .SP per ispezione."""
    with open(path, 'rb') as f:
        raw = f.read()

    if raw[:4] != b'PEPE':
        print("  Non è un file SP valido.")
        return

    desc = raw[4:44].decode('ascii', errors='replace').rstrip('\x00').strip()
    print(f"  Descrizione header: {repr(desc)}")
    print()

    pos = HEADER_SIZE
    blk_idx = 0
    while pos + 6 <= len(raw):
        block_id  = struct.unpack_from('<h', raw, pos)[0]
        block_len = struct.unpack_from('<i', raw, pos + 2)[0]
        bname     = BLOCK_NAMES.get(block_id, f'id={block_id}')
        print(f"  BLOCK[{blk_idx}] {bname}  ({block_len} byte)")
        pos += 6
        block_data = raw[pos:pos + block_len]
        pos += block_len
        blk_idx += 1

        # Prova a parsificare come sequenza di member-block
        mpos = 0
        while mpos + 8 <= len(block_data):
            m_id  = struct.unpack_from('<h', block_data, mpos)[0]
            m_len = struct.unpack_from('<i', block_data, mpos + 2)[0]
            if m_len < 2 or mpos + 6 + m_len > len(block_data):
                break
            typecode = struct.unpack_from('<h', block_data, mpos + 6)[0]
            m_data   = block_data[mpos + 8 : mpos + 6 + m_len]
            mpos    += 6 + m_len

            mname = MEMBER_NAMES.get(m_id, f'id={m_id}')
            text  = _try_decode(m_data)
            extra = f'  → {text}' if text else f'  hex: {m_data[:12].hex(" ")}'
            print(f'    MEMBER {mname:<30s}  typecode={typecode:6d}  {len(m_data):5d}B{extra}')
            if m_id == -29825:  # DataSetHistoryRecord
                decode_history(m_data)

    print()


# ---------------------------------------------------------------------------
def main():
    paths = sys.argv[1:]

    if not paths:
        root = tk.Tk()
        root.withdraw()
        selected = filedialog.askopenfilenames(
            title="Seleziona file .SP PerkinElmer",
            filetypes=[("PerkinElmer FTIR", "*.sp *.SP"), ("Tutti i file", "*.*")]
        )
        root.destroy()
        paths = list(selected)

    if not paths:
        print("Nessun file selezionato. Uscita.")
        return

    fig, ax = plt.subplots(figsize=(10, 5))
    ax.set_xlabel("Wavenumber (cm⁻¹)")
    ax.set_ylabel("Absorbance / Transmittance")
    ax.set_title("Test parser nativo .SP — PerkinElmer FTIR")
    ax.grid(True, linestyle=':', alpha=0.6)

    ok = 0
    for path in paths:
        name = os.path.splitext(os.path.basename(path))[0]
        print(f"\n{'='*60}")
        print(f"File: {path}")
        print(f"{'='*60}")
        dump_sp(path)
        try:
            x, y, desc, date = parse_sp(path)
            print(f"  Descrizione: '{desc}'")
            print(f"  Data misura: {date or 'N/A'}")
            print(f"  X range: {x[0]:.2f} – {x[-1]:.2f}  ({len(x)} punti)")
            print(f"  Y range: {y.min():.4f} – {y.max():.4f}")
            ax.plot(x, y, label=f"{name} ({date or 'N/A'})")
            ok += 1
        except Exception as e:
            print(f"  ERRORE: {e}")

    if ok == 0:
        print("\nNessuno spettro caricato con successo.")
        return

    ax.legend(fontsize=8)
    plt.tight_layout()
    plt.show()


if __name__ == '__main__':
    main()
