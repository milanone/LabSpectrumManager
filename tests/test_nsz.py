"""Horiba SZ-100 `.nsz` (dynamic light scattering) result: reader and metadata attached to a spectrum.

The synthetic files reuse the minimal OLE2 writer of test_vspd_reader; the real-file test compares with the CSV export
of the same measurements and is skipped when the sample is not present.

    py -m unittest discover -s tests -v
"""
import importlib.util
import os
import struct
import tempfile
import tkinter as tk
import unittest
from unittest import mock

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))


def load(name, rel):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, rel))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


lsm = load("lsm_nsz", os.path.join("..", "LabSpectrumManager.pyw"))
build_ole = load("vspd_tests", "test_vspd_reader.py").build_ole

DLS_DIR = os.path.join(HERE, "..", "spettri di esempio", "20261118 lipo vuoti varie size DLS")


def nsz_bytes(values, ints=None, extra=b''):
    """Streams in the layout of a real file: 'Key' + 0 byte + float32 (int32 for the `ints`), split over two streams."""
    def rec(k, v):
        return k.encode() + b'\0' + struct.pack('<f', v) + b'\0\3\x0c\0'

    def rec_int(k, v):
        return k.encode() + b'\0' + struct.pack('<i', v) + b'\0\3\x0c\0'

    calc = b''.join(rec(k, v) for k, v in values.items() if k.startswith(('Cum', 'Calc')))
    meas = (b''.join(rec(k, v) for k, v in values.items() if k.startswith('Meas'))
            + b''.join(rec_int(k, v) for k, v in (ints or {}).items()) + extra)
    tree = {'Base1': calc, 'Base34': meas}
    return build_ole(tree)


GOOD = {'Cum_fMean': 63.0, 'Cum_fPi': 0.1, 'CalcMean[0]': 67.0, 'CalcSD[0]': 17.0, 'CalcMode': 63.5,
        'MeasAngle': 90.0, 'MeasHolderTemp': 25.0, 'MeasSolventVisco': 0.895, 'MeasSolventRef': 1.333,
        'MeasSampleRef': 1.6}


class TestReadNsz(unittest.TestCase):
    def test_values_and_the_radius_to_diameter_conversion(self):
        r = lsm.read_nsz(nsz_bytes(GOOD))
        self.assertAlmostEqual(r['z_average'], 126.0, places=3)        # stored as radius
        self.assertAlmostEqual(r['mean'], 134.0, places=3)
        self.assertAlmostEqual(r['sd'], 34.0, places=3)
        self.assertAlmostEqual(r['mode'], 127.0, places=3)
        self.assertAlmostEqual(r['pi'], 0.1, places=5)
        self.assertAlmostEqual(r['solvent_index'], 1.333, places=5)
        self.assertAlmostEqual(r['particle_index'], 1.6, places=5)
        self.assertEqual(r['angle'], 90.0)

    def test_transmission_before_and_after_are_read_as_integers(self):
        r = lsm.read_nsz(nsz_bytes(GOOD, {'MeasBeforeLightVal': 27454, 'MeasAfterLightVal': 27462}))
        self.assertEqual((r['transmission_before'], r['transmission_after']), (27454.0, 27462.0))

    def test_transmission_out_of_range_is_dropped(self):
        r = lsm.read_nsz(nsz_bytes(GOOD, {'MeasBeforeLightVal': -5, 'MeasAfterLightVal': 27462}))
        self.assertNotIn('transmission_before', r)
        self.assertIn('transmission_after', r)

    def test_missing_optional_keys_are_left_out(self):
        r = lsm.read_nsz(nsz_bytes({'Cum_fMean': 50.0, 'Cum_fPi': 0.2}))
        self.assertEqual(set(r), {'z_average', 'pi'})

    def test_out_of_range_values_are_dropped(self):
        r = lsm.read_nsz(nsz_bytes(dict(GOOD, MeasSolventRef=17.0, MeasAngle=0.0)))
        self.assertNotIn('solvent_index', r)
        self.assertNotIn('angle', r)
        self.assertIn('z_average', r)

    def test_no_z_average_means_not_a_dls_result(self):
        with self.assertRaises(ValueError):
            lsm.read_nsz(nsz_bytes({'MeasAngle': 90.0}))

    def test_not_ole_is_refused(self):
        with self.assertRaises(ValueError):
            lsm.read_nsz(b'plain text' * 40)

    @unittest.skipUnless(os.path.isdir(DLS_DIR), "real DLS sample not present")
    def test_real_files_match_the_csv_export(self):
        # (file, Z-average, PI, mean, S.D., mode) as in the instrument's own CSV of the same measurement
        expected = (("lipo 100 nm 1a10 in PBS_0829.nsz", 125.6, 0.096, 134.2, 34.6, 126.2, 27454, 27462),
                    ("lipo 200 nm 1a40 in PBS_0828.nsz", 177.1, 0.185, 187.0, 42.5, 181.3, 26578, 26534),
                    ("lipo 450 nm 1a40 in PBS_0826.nsz", 299.0, 0.383, 412.8, 96.6, 380.4, 12372, 12529))
        for fn, z, pi, mean, sd, mode, t0, t1 in expected:
            with open(os.path.join(DLS_DIR, fn), 'rb') as f:
                r = lsm.read_nsz(f.read())
            self.assertAlmostEqual(r['z_average'], z, delta=0.06, msg=fn)
            self.assertAlmostEqual(r['pi'], pi, delta=0.0006, msg=fn)
            self.assertAlmostEqual(r['mean'], mean, delta=0.06, msg=fn)
            self.assertAlmostEqual(r['sd'], sd, delta=0.06, msg=fn)
            self.assertAlmostEqual(r['mode'], mode, delta=0.06, msg=fn)
            # transmitted intensity before/after: the values of the CSV export and of the instrument window (photos)
            self.assertEqual((r['transmission_before'], r['transmission_after']), (t0, t1), msg=fn)
            self.assertEqual((round(r['solvent_index'], 3), round(r['particle_index'], 2), r['angle']), (1.333, 1.6, 90.0))


class TestAttach(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        try:
            cls.root = lsm.TkinterDnD.Tk() if lsm.HAS_DND else tk.Tk()
        except tk.TclError as e:
            raise unittest.SkipTest("Tk not available: %s" % e)
        cls.root.withdraw()

    @classmethod
    def tearDownClass(cls):
        cls.root.destroy()

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        with mock.patch.object(lsm.messagebox, "showwarning"):
            self.app = lsm.LabSpectrumManager(self.root)
        wl = np.arange(350, 801)
        path = os.path.join(self.tmp.name, "lip.csv")
        with open(path, "w", newline="") as f:
            f.write("Wavelength [nm],Sample X1\n")
            for w in wl:
                f.write(f"{w},{0.07 * (w / 1000.0) ** -2.6:.6f}\n")
        self.app.processa_file(path)
        self.app.aggiorna_vista()
        self.name = next(iter(self.app.spectra))
        self.nsz = os.path.join(self.tmp.name, "r.nsz")
        with open(self.nsz, 'wb') as f:
            f.write(nsz_bytes(GOOD, {'MeasBeforeLightVal': 27454, 'MeasAfterLightVal': 27462}))

    def test_attach_adds_the_dls_rows_to_the_metadata(self):
        self.assertTrue(self.app._allega_dls_file(self.name, self.nsz))
        info = self.app.spectra[self.name]['info']
        self.assertEqual(info['DLS Z-average (nm)'], '126.0')
        self.assertEqual(info['DLS PI'], '0.100')
        self.assertEqual(info['DLS file'], 'r.nsz')
        self.assertEqual((info['DLS transmission before'], info['DLS transmission after']), ('27454', '27462'))
        self.assertIn('DLS Z-average (nm)', self.app.param_box.get('1.0', tk.END))
        self.assertTrue(self.app._dirty)

    def test_no_dls_means_no_dls_rows(self):
        info = self.app.spectra[self.name]['info']
        self.assertFalse([k for k in info if k.startswith('DLS')])
        self.assertNotIn('DLS', self.app.param_box.get('1.0', tk.END))

    def test_a_second_attachment_replaces_the_first(self):
        self.app._allega_dls_file(self.name, self.nsz)
        other = os.path.join(self.tmp.name, "r2.nsz")
        with open(other, 'wb') as f:
            f.write(nsz_bytes({'Cum_fMean': 90.0, 'Cum_fPi': 0.2}))      # no viscosity etc.
        self.app._allega_dls_file(self.name, other)
        info = self.app.spectra[self.name]['info']
        self.assertEqual(info['DLS Z-average (nm)'], '180.0')
        self.assertNotIn('DLS viscosity (mPa.s)', info)                   # nothing left from the first one

    def test_a_bad_file_changes_nothing_and_reports(self):
        bad = os.path.join(self.tmp.name, "bad.nsz")
        with open(bad, 'wb') as f:
            f.write(b'not a dls file' * 30)
        with mock.patch.object(lsm.messagebox, "showerror") as err:
            self.assertFalse(self.app._allega_dls_file(self.name, bad))
        err.assert_called_once()
        self.assertFalse([k for k in self.app.spectra[self.name]['info'] if k.startswith('DLS')])

    @unittest.skipUnless(lsm.__dict__.get('HAS_SCIPY', True), "n/a")
    def test_simulator_starts_from_the_dls_size_and_width_or_from_the_defaults(self):
        try:
            import scipy.special  # noqa: F401
        except ImportError:
            self.skipTest("scipy not installed")
        before = set(self.root.winfo_children())
        self.app.open_simulator()
        win = [w for w in self.root.winfo_children() if w not in before and isinstance(w, tk.Toplevel)][0]
        self.addCleanup(win.destroy)
        self.assertEqual((win._sim_vars['d'].get(), win._sim_vars['sg'].get()), ('150', '0.2'))      # no DLS: defaults
        win.destroy()
        self.app._allega_dls_file(self.name, self.nsz)
        before = set(self.root.winfo_children())
        self.app.open_simulator()
        win = [w for w in self.root.winfo_children() if w not in before and isinstance(w, tk.Toplevel)][0]
        self.addCleanup(win.destroy)
        self.assertEqual(win._sim_vars['d'].get(), '126')                                            # Z-average
        self.assertEqual(win._sim_vars['sg'].get(), '0.32')                                          # sqrt(PI 0.1)


if __name__ == "__main__":
    unittest.main()
