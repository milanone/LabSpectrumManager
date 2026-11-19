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


def nsz_bytes(values, ints=None, extra=b'', correlogram=None):
    """Streams in the layout of a real file: 'Key' + 0 byte + float32 (int32 for the `ints`), split over two streams.
    `correlogram` = (delays in s, normalised ACF): stored as the 'Delay' and 'OriginalACF' arrays of a third stream, each
    preceded by its byte length, the ACF followed by three junk values as in the real files."""
    def rec(k, v):
        return k.encode() + b'\0' + struct.pack('<f', v) + b'\0\3\x0c\0'

    def rec_int(k, v):
        return k.encode() + b'\0' + struct.pack('<i', v) + b'\0\3\x0c\0'

    def rec_arr(k, a):
        return b'\3\6\0' + struct.pack('<I', 4 * len(a)) + k.encode() + b'\0' + struct.pack('<%df' % len(a), *a)
    calc = b''.join(rec(k, v) for k, v in values.items() if k.startswith(('Cum', 'Calc')))
    meas = (b''.join(rec(k, v) for k, v in values.items() if k.startswith('Meas'))
            + b''.join(rec_int(k, v) for k, v in (ints or {}).items()) + extra)
    tree = {'Base1': calc, 'Base34': meas}
    if correlogram is not None:
        tau, acf = correlogram
        tree['Base2'] = rec_arr('Delay', tau) + rec_arr('OriginalACF', list(acf) + [1e-4, 0.0, 1.0e9])
    return build_ole(tree)


def synthetic_correlogram(populations, q2=None, T=25.0, visc=0.895, n=1.333, angle=90.0):
    """g1 of an intensity-weighted mixture of narrow log-normal peaks [(median nm, sigma, fraction), ...]."""
    q2 = lsm.dls_q2(angle, n) if q2 is None else q2
    tau = np.geomspace(1e-7, 0.11, 43)
    d = np.geomspace(5.0, 5000.0, 400)
    p = np.zeros_like(d)
    for med, sg, frac in populations:
        g = np.exp(-0.5 * (np.log(d / med) / sg) ** 2)
        p += frac * g / g.sum()
    g1 = np.exp(-np.outer(tau, lsm._dls_gamma(d, q2, T, visc))) @ p
    return tau, g1 ** 2


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


try:
    import scipy.optimize  # noqa: F401
    HAS_SCIPY = True
except ImportError:
    HAS_SCIPY = False


class TestCorrelogram(unittest.TestCase):
    def rec(self, tau, acf):
        return {'tau': tau, 'acf': acf, 'angle': 90.0, 'n': 1.333, 'visc': 0.895, 'T': 25.0}

    def test_the_arrays_are_read_from_the_file(self):
        tau, acf = synthetic_correlogram([(100.0, 0.05, 1.0)])
        r = lsm.read_nsz(nsz_bytes(GOOD, correlogram=(tau, acf)))
        np.testing.assert_allclose(r['tau'], tau, rtol=1e-6)
        np.testing.assert_allclose(r['acf'], acf, atol=1e-6)          # the three junk values are cut off
        self.assertEqual(len(r['acf']), 43)

    def test_a_file_without_a_correlogram_has_no_arrays(self):
        r = lsm.read_nsz(nsz_bytes(GOOD))
        self.assertNotIn('tau', r)

    def test_cumulants_recover_the_diameter(self):
        tau, acf = synthetic_correlogram([(100.0, 0.05, 1.0)])
        z, pi = lsm.dls_cumulants(tau, acf, lsm.dls_q2(90.0, 1.333), 25.0, 0.895)
        self.assertAlmostEqual(z, 100.0, delta=3.0)
        self.assertLess(pi, 0.05)

    @unittest.skipUnless(HAS_SCIPY, "scipy not installed")
    def test_regularised_inversion_finds_the_mode(self):
        tau, acf = synthetic_correlogram([(100.0, 0.08, 1.0)])
        r = lsm.dls_analyse(self.rec(tau, acf), lsm.DLS_MODELS[0])
        self.assertAlmostEqual(r['mode'], 100.0, delta=15.0)

    @unittest.skipUnless(HAS_SCIPY, "scipy not installed")
    def test_one_peak_fit_recovers_median_and_width(self):
        tau, acf = synthetic_correlogram([(150.0, 0.2, 1.0)])
        r = lsm.dls_analyse(self.rec(tau, acf), lsm.DLS_MODELS[1])
        pk = r['peaks'][0]
        self.assertAlmostEqual(pk['median'], 150.0, delta=10.0)
        self.assertAlmostEqual(pk['sigma'], 0.2, delta=0.06)

    @unittest.skipUnless(HAS_SCIPY, "scipy not installed")
    def test_two_peaks_are_resolved_and_preferred_by_the_bic_when_they_are_far_apart(self):
        tau, acf = synthetic_correlogram([(40.0, 0.08, 0.3), (400.0, 0.08, 0.7)])
        rec = self.rec(tau, acf)
        r2, r1 = lsm.dls_analyse(rec, lsm.DLS_MODELS[2]), lsm.dls_analyse(rec, lsm.DLS_MODELS[1])
        small, large = r2['peaks']
        self.assertAlmostEqual(small['median'], 40.0, delta=10.0)
        self.assertAlmostEqual(large['median'], 400.0, delta=60.0)
        self.assertAlmostEqual(small['fraction'], 0.3, delta=0.1)
        self.assertLess(r2['bic'], r1['bic'] - 6)

    @unittest.skipUnless(HAS_SCIPY, "scipy not installed")
    def test_a_single_peak_is_not_split_by_the_bic(self):
        tau, acf = synthetic_correlogram([(200.0, 0.15, 1.0)])
        rec = self.rec(tau, acf)
        r2, r1 = lsm.dls_analyse(rec, lsm.DLS_MODELS[2]), lsm.dls_analyse(rec, lsm.DLS_MODELS[1])
        self.assertGreater(r2['bic'], r1['bic'] - 6)           # two peaks bring no significant improvement

    @unittest.skipUnless(os.path.isdir(DLS_DIR) and HAS_SCIPY, "real DLS sample or scipy not present")
    def test_real_files_reproduce_the_instrument_z_average_and_mode(self):
        for fn, z, mode in (("lipo 100 nm 1a10 in PBS_0829.nsz", 125.6, 134.2), ("lipo 200 nm 1a40 in PBS_0828.nsz", 177.1, 193.5)):
            with open(os.path.join(DLS_DIR, fn), 'rb') as f:
                r = lsm.read_nsz(f.read())
            rec = {'tau': r['tau'], 'acf': r['acf'], 'angle': r['angle'], 'n': r['solvent_index'], 'visc': r['viscosity'],
                   'T': r['temperature']}
            a = lsm.dls_analyse(rec, lsm.DLS_MODELS[0])
            self.assertAlmostEqual(a['z'], z, delta=0.03 * z, msg=fn)           # the cumulant Z-average, within 3%
            self.assertAlmostEqual(a['mode'], mode, delta=0.08 * mode, msg=fn)  # the mode of the instrument's table, 8%


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
        self.app.apri_correzione_scattering()
        before = set(self.root.winfo_children())
        self.app._sc_simulatore()
        win = [w for w in self.root.winfo_children() if w not in before and isinstance(w, tk.Toplevel)][0]
        self.addCleanup(win.destroy)
        self.assertEqual((win._sim_vars['d'].get(), win._sim_vars['sg'].get()), ('150', '0.2'))      # no DLS: defaults
        win.destroy()
        self.app._sc_annulla()
        self.app._allega_dls_file(self.name, self.nsz)
        self.app.apri_correzione_scattering()
        before = set(self.root.winfo_children())
        self.app._sc_simulatore()
        win = [w for w in self.root.winfo_children() if w not in before and isinstance(w, tk.Toplevel)][0]
        self.addCleanup(win.destroy)
        self.assertEqual(win._sim_vars['d'].get(), '126')                                            # Z-average
        self.assertEqual(win._sim_vars['sg'].get(), '0.32')                                          # sqrt(PI 0.1)


@unittest.skipUnless(HAS_SCIPY, "scipy not installed")
class TestCorrelogramWindow(TestAttach):
    def attach_with_correlogram(self):
        tau, acf = synthetic_correlogram([(60.0, 0.08, 0.3), (300.0, 0.08, 0.7)])
        path = os.path.join(self.tmp.name, "c.nsz")
        with open(path, 'wb') as f:
            f.write(nsz_bytes(GOOD, correlogram=(tau, acf)))
        self.assertTrue(self.app._allega_dls_file(self.name, path))

    def open_window(self):
        before = set(self.root.winfo_children())
        self.app._dls_finestra(self.name)
        win = [w for w in self.root.winfo_children() if w not in before and isinstance(w, tk.Toplevel)][0]
        self.addCleanup(win.destroy)
        return win

    @staticmethod
    def label_text(win):
        return [w.cget("text") for w in win.winfo_children() if isinstance(w, tk.Label)][-1]

    def test_the_attachment_stores_the_correlogram(self):
        self.attach_with_correlogram()
        rec = self.app.spectra[self.name]['dls']
        self.assertEqual(len(rec['tau']), 43)
        self.assertEqual((rec['angle'], rec['T']), (90.0, 25.0))

    def test_a_new_attachment_without_a_correlogram_drops_the_old_one(self):
        self.attach_with_correlogram()
        self.assertTrue(self.app._allega_dls_file(self.name, self.nsz))        # the plain file of setUp
        self.assertNotIn('dls', self.app.spectra[self.name])

    def test_window_without_a_correlogram_warns_and_does_not_open(self):
        with mock.patch.object(lsm.messagebox, "showwarning") as warn:
            before = set(self.root.winfo_children())
            self.app._dls_finestra(self.name)
        warn.assert_called_once()
        self.assertEqual(set(self.root.winfo_children()), before)

    def test_window_reports_each_model(self):
        self.attach_with_correlogram()
        win = self.open_window()
        self.assertIn("Distribution: mean", self.label_text(win))                 # regularised, the first model
        win._dls_mode.set(lsm.DLS_MODELS[2])
        win._dls_compute()
        text = self.label_text(win)
        self.assertIn("Peak 1:", text)
        self.assertIn("Peak 2:", text)
        self.assertIn("two peaks are supported", text)                            # the BIC verdict
        win._dls_mode.set(lsm.DLS_MODELS[1])
        win._dls_compute()
        self.assertIn("Peak 1:", self.label_text(win))
        self.assertNotIn("Peak 2:", self.label_text(win))


if __name__ == "__main__":
    unittest.main()
