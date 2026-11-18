"""Tests for the Shimadzu UV-1900 .vspd reader (and the shared UVProbe .SPC path).

The files are synthetic: a minimal OLE2 (compound file) writer builds the same stream layout
as a real LabSolutions UV-Vis file, so no real data and no Tk window are needed. The test
against the real file is skipped when the sample is not present.

    py -m unittest discover -s tests -v
"""
import importlib.util
import os
import struct
import tempfile
import unittest

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location(
    "lsm", os.path.join(HERE, "..", "LabSpectrumManager.pyw"))
lsm = importlib.util.module_from_spec(spec)
spec.loader.exec_module(lsm)

REAL_FILE = os.path.join(HERE, "..", "spettri di esempio", "lipo bupropione 1a10 in PBS.vspd")

SIG = b'\xd0\xcf\x11\xe0\xa1\xb1\x1a\xe1'
NOSTREAM, ENDOFCHAIN, FATSECT = 0xFFFFFFFF, 0xFFFFFFFE, 0xFFFFFFFD
SEC, MINI, CUTOFF = 512, 64, 4096
DATASET = 'DataStorage1/DataSetGroup/DataSet1/'


def build_ole(tree):
    """Minimal OLE2 v3 writer. `tree` is a nested dict {name: bytes | dict}. Streams shorter than
    4096 bytes go to the mini stream, the others to regular sectors, as in real files (so both
    reader paths are exercised). Siblings are chained through 'right' only, which the reader's
    tree walk accepts."""
    entries = [None]

    def add(nodes):
        idxs = []
        for name, val in nodes.items():
            i = len(entries)
            entries.append(None)
            if isinstance(val, dict):
                entries[i] = dict(name=name, kind=1, child=add(val), data=b'')
            else:
                entries[i] = dict(name=name, kind=2, child=NOSTREAM, data=val)
            entries[i].update(left=NOSTREAM, right=NOSTREAM, start=ENDOFCHAIN)
            idxs.append(i)
        for a, b in zip(idxs, idxs[1:]):
            entries[a]['right'] = b
        return idxs[0] if idxs else NOSTREAM

    root_child = add(tree)
    mini, minifat = bytearray(), []
    for e in entries[1:]:
        if e['kind'] == 2 and len(e['data']) < CUTOFF:
            n = -(-len(e['data']) // MINI)
            e['start'] = len(minifat)
            for k in range(n):
                minifat.append(len(minifat) + 1 if k < n - 1 else ENDOFCHAIN)
            mini += e['data'].ljust(n * MINI, b'\0')
    entries[0] = dict(name='Root Entry', kind=5, child=root_child, data=bytes(mini),
                      left=NOSTREAM, right=NOSTREAM, start=ENDOFCHAIN)

    def sectors(n_bytes):
        return -(-n_bytes // SEC)

    n_dir = sectors(-(-len(entries) // 4) * 4 * 128)
    n_minifat = sectors(len(minifat) * 4)
    blobs = [('dir', None), ('minifat', None), ('root', entries[0])]
    blobs += [('big', e) for e in entries[1:] if e['kind'] == 2 and len(e['data']) >= CUTOFF]
    sizes = {'dir': n_dir, 'minifat': n_minifat, 'root': sectors(len(mini))}
    total = sum(sizes[k] if k in sizes else sectors(len(e['data'])) for k, e in blobs)
    n_fat = 1
    while n_fat * 128 < total + n_fat:
        n_fat += 1

    fat = [NOSTREAM] * (n_fat * 128)
    for i in range(n_fat):
        fat[i] = FATSECT
    nxt = n_fat
    placed = {}
    for kind, e in blobs:
        n = sizes[kind] if kind in sizes else sectors(len(e['data']))
        if n == 0:
            continue
        for k in range(n):
            fat[nxt + k] = nxt + k + 1 if k < n - 1 else ENDOFCHAIN
        placed[(kind, id(e))] = nxt
        nxt += n
    if n_minifat == 0:
        minifat_start = ENDOFCHAIN
    else:
        minifat_start = placed[('minifat', id(None))]

    dir_start = placed[('dir', id(None))]
    root_start = placed[('root', id(entries[0]))]
    entries[0]['start'] = root_start
    for e in entries[1:]:
        if e['kind'] == 2 and len(e['data']) >= CUTOFF:
            e['start'] = placed[('big', id(e))]

    def dir_entry(e):
        if e is None:
            return struct.pack('<64sHBB3I16sI16sIQ', b'', 0, 0, 1, NOSTREAM, NOSTREAM, NOSTREAM,
                               b'', 0, b'', 0, 0)
        raw = (e['name'] + '\0').encode('utf-16le')
        size = len(e['data']) if e['kind'] != 1 else 0
        return struct.pack('<64sHBB3I16sI16sIQ', raw, len(raw), e['kind'], 1,
                           e['left'], e['right'], e['child'], b'', 0, b'', e['start'], size)

    directory = b''.join(dir_entry(e) for e in entries)
    directory = directory.ljust(n_dir * SEC, b'\0')
    body = {}
    body['dir'] = directory
    body['minifat'] = struct.pack('<%dI' % len(minifat), *minifat).ljust(n_minifat * SEC, b'\0')
    body['root'] = bytes(mini).ljust(sizes['root'] * SEC, b'\0')

    out = bytearray(SEC * (1 + nxt))
    for kind, e in blobs:
        data = body[kind] if kind in body else e['data']
        data = data.ljust(sectors(len(data)) * SEC, b'\0')
        if not data:
            continue
        start = placed[(kind, id(e) if kind != 'dir' and kind != 'minifat' else id(None))]
        out[(start + 1) * SEC:(start + 1) * SEC + len(data)] = data
    fat_bytes = struct.pack('<%dI' % len(fat), *fat)
    out[SEC:SEC + len(fat_bytes)] = fat_bytes

    difat = list(range(n_fat)) + [NOSTREAM] * (109 - n_fat)
    header = (SIG + b'\0' * 16 + struct.pack('<HHHHH6x', 0x3E, 3, 0xFFFE, 9, 6)
              + struct.pack('<IIIIIIIII', 0, n_fat, dir_start, 0, CUTOFF, minifat_start,
                            n_minifat, ENDOFCHAIN, 0)
              + struct.pack('<109I', *difat))
    out[:len(header)] = header
    return bytes(out)


def vspd_bytes(x, y, date=45790.5, kind='RawData', instrument='UV-1900 Series'):
    """Same stream names and encodings as the real file: float64 X/Y, a date after a
    length-prefixed string in DataSetHeaderInfo, a length-prefixed instrument name."""
    hdr = bytes([len(kind)]) + kind.encode() + b'\0' * 4 + struct.pack('<d', date) + b'\0' * 40
    ins = struct.pack('<I', 17) + bytes([len(instrument)]) + instrument.encode()
    return build_ole({'DataStorage1': {'DataSetGroup': {'DataSet1': {
        'DataSetHeaderInfo': hdr,
        'DataSpectrumStorage': {'Data': {
            'X Data.1': np.asarray(x, '<f8').tobytes(),
            'Y Data.1': np.asarray(y, '<f8').tobytes()}},
        'MethodStorage': {'InstrumentNameStream': ins},
    }}}})


class VspdReader(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.app = lsm.LabSpectrumManager.__new__(lsm.LabSpectrumManager)
        self.app.spectra = {}

    def write(self, name, data):
        p = os.path.join(self.tmp.name, name)
        with open(p, 'wb') as f:
            f.write(data)
        return p

    def test_synthetic_file_values_and_order(self):
        x = np.arange(600, 199, -1.0)                       # scan order: descending, 401 points
        y = 0.8 * (x / 1000.0) ** -2.5 * 0.1 + 0.02
        path = self.write('lipo test.vspd', vspd_bytes(x, y))
        self.app.processa_file(path)
        self.assertEqual(list(self.app.spectra), ['lipo test'])
        s = self.app.spectra['lipo test']
        df = s['df']
        self.assertEqual(len(df), 401)
        self.assertTrue(df.index.is_monotonic_increasing)   # re-sorted ascending
        self.assertEqual((df.index.min(), df.index.max()), (200.0, 600.0))
        np.testing.assert_allclose(df['lipo test'].to_numpy(), y[::-1], rtol=5e-4)  # 4 significant digits

    def test_metadata(self):
        x = np.arange(600, 199, -1.0)
        path = self.write('m.vspd', vspd_bytes(x, np.full(x.size, 0.5)))
        self.app.processa_file(path)
        info = self.app.spectra['m']['info']
        self.assertEqual(info['Type'], 'UV-Vis')
        self.assertEqual(info['Date'], '13/05/25')            # 45790.5 = 2025-05-13 12:00
        self.assertEqual(info['Instrument'], 'Shimadzu UV-1900 Series')
        self.assertEqual(info['Scan Range'], '200.0 - 600.0 nm')

    def test_date_offset_follows_the_string_length(self):
        """'RawData' (.vspd) and 'Manipulated' (.SPC) put the date at different offsets."""
        x = np.arange(600, 199, -1.0)
        path = self.write('n.vspd', vspd_bytes(x, np.ones(x.size), kind='Manipulated'))
        self.app.processa_file(path)
        self.assertEqual(self.app.spectra['n']['info']['Date'], '13/05/25')

    def test_not_ole_is_refused(self):
        path = self.write('bad.vspd', b'not an ole file' * 50)
        with self.assertRaises(ValueError):
            self.app.leggi_vspd(path)

    @unittest.skipUnless(os.path.exists(REAL_FILE), "real .vspd sample not present")
    def test_real_file(self):
        """Values read independently with olefile while decoding the format (not from a
        LabSolutions export: none was available for this file)."""
        self.app.leggi_vspd(REAL_FILE)
        s = self.app.spectra['lipo bupropione 1a10 in PBS']
        df = s['df']
        col = df.iloc[:, 0]
        self.assertEqual(len(df), 401)
        self.assertEqual((df.index.min(), df.index.max()), (200.0, 600.0))
        self.assertAlmostEqual(col.loc[600.0], 0.0701, places=3)
        self.assertAlmostEqual(col.loc[400.0], 0.2184, places=3)
        self.assertAlmostEqual(col.loc[200.0], 1.9591, places=3)
        self.assertEqual(s['info']['Date'], '13/05/25')


if __name__ == '__main__':
    unittest.main()
