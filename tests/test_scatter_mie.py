"""Pure power-law fit (offset fixed at 0) and the theoretical P(d) curve of vesicles (Mie theory).

    python -m unittest discover -s tests -v
"""
import importlib.util
import os
import tempfile
import tkinter as tk
import unittest
from unittest import mock

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("lsm_mie", os.path.join(HERE, "..", "LabSpectrumManager.pyw"))
lsm = importlib.util.module_from_spec(spec)
spec.loader.exec_module(lsm)

try:
    import scipy.special  # noqa: F401
    HAS_SCIPY = True
except ImportError:
    HAS_SCIPY = False


class TestFitPowerLaw(unittest.TestCase):
    def setUp(self):
        self.x = np.arange(350.0, 801.0, 5.0)

    def test_recovers_exponent_and_scale(self):
        y = 0.08 * (self.x / 1000.0) ** -2.7
        vs, p, m, rms = lsm.fit_power_law(self.x, y)
        self.assertAlmostEqual(p, 2.7, places=2)
        self.assertAlmostEqual(vs, 0.08, places=4)
        self.assertLess(rms, 1e-5)

    def test_fixed_exponent_and_fixed_slope(self):
        y = 0.05 * (self.x / 1000.0) ** -3.0 - 2e-5 * self.x
        vs, p, m, rms = lsm.fit_power_law(self.x, y, P=3.0, m=2e-5)
        self.assertAlmostEqual(vs, 0.05, places=5)
        self.assertEqual(p, 3.0)

    def test_free_slope(self):
        y = 0.05 * (self.x / 1000.0) ** -3.0 - 2e-5 * self.x
        vs, p, m, rms = lsm.fit_power_law(self.x, y, m=None)
        self.assertAlmostEqual(p, 3.0, places=1)
        self.assertAlmostEqual(m, 2e-5, delta=3e-6)

    def test_offset_is_not_absorbed(self):
        """With a real offset the pure power law must NOT hide it: the residual shows it."""
        y = 0.05 * (self.x / 1000.0) ** -3.0 + 0.02
        _, _, _, rms = lsm.fit_power_law(self.x, y)
        self.assertGreater(rms, 1e-3)


@unittest.skipUnless(HAS_SCIPY, "scipy not installed")
class TestVesicleTheory(unittest.TestCase):
    lam = np.arange(350.0, 801.0, 15.0)

    def test_rayleigh_limit_gives_P_of_four(self):
        P = lsm.vesicle_P_curve([20.0], self.lam, sigma=0.0, acc_deg=0.0)[0]
        self.assertAlmostEqual(P, 4.0, delta=0.05)

    def test_P_falls_with_size_in_the_useful_range(self):
        P = lsm.vesicle_P_curve([40.0, 80.0, 120.0], self.lam, sigma=0.2, acc_deg=10.0)
        self.assertTrue(P[0] > P[1] > P[2])
        self.assertGreater(P[0], 3.4)
        self.assertLess(P[2], 3.0)

    def test_collection_angle_reduces_the_extinction_of_large_particles(self):
        full = lsm.vesicle_turbidity(self.lam, 400.0, acc_deg=0.0)
        part = lsm.vesicle_turbidity(self.lam, 400.0, acc_deg=10.0)
        self.assertTrue(np.all(part < full))
        self.assertTrue(np.all(part > 0))

    def test_coated_sphere_equals_the_thin_shell_limit_of_a_solid_one(self):
        """A shell as thick as the radius is a solid sphere: compare with the textbook Mie sum."""
        lam = np.array([500.0])
        d, dn = 300.0, 0.05
        got = lsm.vesicle_turbidity(lam, d, t=d / 2, dn=dn, acc_deg=0.0)[0]
        from scipy.special import spherical_jn, spherical_yn
        nm = lsm._water_index(lam)[0]
        m = (nm + dn) / nm
        x = 2 * np.pi * nm / 500.0 * d / 2
        s = 0.0
        for n in range(1, int(x + 4 * x ** (1 / 3) + 2) + 1):
            def ric(z):
                j, y = spherical_jn(n, z), spherical_yn(n, z)
                return z * j, j + z * spherical_jn(n, z, derivative=True), \
                    z * (j - 1j * y), (j - 1j * y) + z * (spherical_jn(n, z, derivative=True)
                                                           - 1j * spherical_yn(n, z, derivative=True))
            pm, dm, _, _ = ric(m * x)
            p, dp, xi, dxi = ric(x)
            a = (m * pm * dp - p * dm) / (m * pm * dxi - xi * dm)
            b = (pm * dp - m * p * dm) / (pm * dxi - m * xi * dm)
            s += (2 * n + 1) * (abs(a) ** 2 + abs(b) ** 2)
        ref = 2 / x ** 2 * s * np.pi * (d / 2) ** 2
        # the coated-sphere code keeps a 1 nm hole at the centre; relative difference stays tiny
        self.assertAlmostEqual(got / ref, 1.0, delta=0.01)


@unittest.skipUnless(HAS_SCIPY, "scipy not installed")
class TestAmplitude(unittest.TestCase):
    lam = np.array([500.0])

    def test_unilamellar_reference_value(self):
        """Value computed while developing the model (100 nm, 10 mg/mL POPC, 0.1 cm cuvette)."""
        a = lsm.vesicle_amplitude(self.lam, 100.0, 0.2, 1, 10.0, 0.1)[0]
        self.assertAlmostEqual(a, 0.0606, delta=0.002)

    def test_proportional_to_concentration_and_path(self):
        a = lsm.vesicle_amplitude(self.lam, 150.0, 0.2, 1, 10.0, 1.0)[0]
        self.assertAlmostEqual(lsm.vesicle_amplitude(self.lam, 150.0, 0.2, 1, 2.0, 1.0)[0] / a, 0.2, places=3)
        self.assertAlmostEqual(lsm.vesicle_amplitude(self.lam, 150.0, 0.2, 1, 10.0, 0.1)[0] / a, 0.1, places=3)

    def test_amplitude_grows_about_linearly_with_the_lamellae(self):
        a1 = lsm.vesicle_amplitude(self.lam, 400.0, 0.3, 1, 0.25, 1.0)[0]
        a7 = lsm.vesicle_amplitude(self.lam, 400.0, 0.3, 7, 0.25, 1.0)[0]
        self.assertTrue(6.0 < a7 / a1 < 9.0)

    def test_amplitude_grows_with_size_at_fixed_lipid(self):
        a = [lsm.vesicle_amplitude(self.lam, d, 0.2, 1, 10.0, 0.1)[0] for d in (60.0, 100.0, 200.0)]
        self.assertTrue(a[0] < a[1] < a[2])


@unittest.skipUnless(HAS_SCIPY, "scipy not installed")
class TestComposition(unittest.TestCase):
    def test_popc_alone_changes_nothing(self):
        mx = lsm.lipid_mixture()
        self.assertEqual((mx['mw'], mx['area'], mx['t'], mx['dn'], mx['peg_film']), (760.08, 0.65, 4.5, 0.12, 0.0))

    def test_85_10_5_mixture(self):
        mx = lsm.lipid_mixture(10, 5)
        self.assertAlmostEqual(mx['mw'], 825.0, delta=0.5)            # 0.85*760.08 + 0.10*386.65 + 0.05*2805.5
        self.assertAlmostEqual(mx['area'], 0.5915, places=3)
        self.assertAlmostEqual(mx['t'], 4.9, places=6)
        self.assertAlmostEqual(mx['peg_film'], 0.234, delta=0.01)     # nm of pure PEG per leaflet
        self.assertGreater(mx['dn'], 0.12 + 0.01)                     # the PEG film adds polarisability

    def test_fractions_above_100_percent_are_refused(self):
        with self.assertRaises(ValueError):
            lsm.lipid_mixture(60, 50)

    def test_no_peg_leaflets_removes_the_film_term(self):
        a = lsm.lipid_mixture(10, 5, peg_leaflets=0)
        b = lsm.lipid_mixture(10, 5, peg_leaflets=2)
        self.assertAlmostEqual(a['dn'], 0.12, places=9)
        self.assertGreater(b['dn'], a['dn'])


class TestPanelOffsetLock(unittest.TestCase):
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
        y = 0.07 * (wl / 1000.0) ** -2.6
        path = os.path.join(self.tmp.name, "lip.csv")
        with open(path, "w", newline="") as f:
            f.write("Wavelength [nm],Sample X1\n")
            for w, v in zip(wl, y):
                f.write(f"{w},{v:.6f}\n")
        self.app.processa_file(path)
        self.app.aggiorna_vista()
        self.app.apri_correzione_scattering()

    def test_auto_with_offset_fixed_returns_the_true_exponent_and_zero_offset(self):
        self.app._sc_OFF_lock.set(True)
        self.app._sc_auto()
        self.assertAlmostEqual(self.app._sc_P['value'].get(), 2.6, places=2)
        self.assertEqual(self.app._sc_OFF['value'].get(), 0.0)
        self.assertEqual(self.app._sc_M['value'].get(), 0.0)           # slope stays fixed at 0 by default
        self.assertIn("offset fixed at 0", self.app._sc_fit_label.cget("text"))

    def test_ticking_the_box_zeroes_the_offset(self):
        self.app._linked_set(self.app._sc_OFF, 0.05)
        self.app._sc_OFF_lock.set(True)
        self.app._sc_OFF_lock.get() and self.app._linked_set(self.app._sc_OFF, 0.0)   # what the checkbox command does
        self.assertEqual(self.app._sc_OFF['value'].get(), 0.0)

    @unittest.skipUnless(HAS_SCIPY, "scipy not installed")
    def test_theory_window_opens_and_reports_a_crossing(self):
        self.app._sc_OFF_lock.set(True)
        self.app._sc_auto()
        before = set(self.root.winfo_children())
        self.app._sc_curva_teorica()
        win = [w for w in self.root.winfo_children() if w not in before and isinstance(w, tk.Toplevel)][0]
        self.addCleanup(win.destroy)
        labels = [w for w in win.winfo_children() if isinstance(w, tk.Label) and "σ" in str(w.cget("text"))]
        self.assertTrue(labels)
        self.assertIn("meets the curve", labels[0].cget("text"))


    def open_simulator(self, x, y):
        """Open the simulator on the spectrum (x, y) as the panel holds it; returns the window and its field variables."""
        self.app._sc_x, self.app._sc_y = np.asarray(x, float), np.asarray(y, float)
        before = set(self.root.winfo_children())
        self.app._sc_simulatore()
        win = [w for w in self.root.winfo_children() if w not in before and isinstance(w, tk.Toplevel)][0]
        self.addCleanup(win.destroy)
        return win

    @unittest.skipUnless(HAS_SCIPY, "scipy not installed")
    def test_simulator_overlays_theory_and_reports_the_agreement(self):
        lam = np.arange(350.0, 801.0, 5.0)
        truth = lsm.vesicle_amplitude(lam, 120.0, 0.2, 1, 10.0 / 40, 1.0)          # 10 mg/mL, 1:40, 1 cm
        win = self.open_simulator(lam, truth)
        v = win._sim_vars
        v['dil'].set('40'); v['d'].set('120')
        win._sim_update()
        text = [w.cget("text") for w in win.winfo_children() if isinstance(w, tk.Label)][-1]
        self.assertIn("rms residual", text)
        rms = float(text.split("rms residual")[1].split()[0])
        self.assertLess(rms, 1e-3 * truth.max())                    # same model, same parameters: no residual

    @unittest.skipUnless(HAS_SCIPY, "scipy not installed")
    def test_simulator_composition_fields_change_the_theory(self):
        lam = np.arange(350.0, 801.0, 10.0)
        truth = lsm.vesicle_amplitude(lam, 150.0, 0.2, 1, 10.0 / 40, 1.0)
        win = self.open_simulator(lam, truth)
        v = win._sim_vars
        v['dil'].set('40'); v['d'].set('150')
        win._sim_update()
        rms0 = float([w.cget("text") for w in win.winfo_children() if isinstance(w, tk.Label)][-1]
                     .split("rms residual")[1].split()[0])
        v['chol'].set('10'); v['peg'].set('5')
        win._sim_update()
        text = [w.cget("text") for w in win.winfo_children() if isinstance(w, tk.Label)][-1]
        self.assertIn("Mixture 85/10/5 mol%", text)
        rms1 = float(text.split("rms residual")[1].split()[0])
        self.assertGreater(rms1, 5 * max(rms0, 1e-6))               # a different bilayer no longer matches POPC-only data

    @unittest.skipUnless(HAS_SCIPY, "scipy not installed")
    def test_every_simulator_field_has_a_help_text_that_can_be_shown(self):
        lam = np.arange(350.0, 801.0, 10.0)
        win = self.open_simulator(lam, lsm.vesicle_amplitude(lam, 150.0, 0.2, 1, 10.0, 1.0))
        self.assertTrue(set(win._sim_vars) <= set(lsm.SIM_TIPS))
        self.assertEqual(len(win._sim_tips), 2 * len(win._sim_vars) + 3)          # label + entry per field, 3 buttons
        tip = next(t for t in win._sim_tips if t.text == lsm.SIM_TIPS['sg'])
        tip._show()
        self.assertIsNotNone(tip._tip)
        self.assertIn("PDI", tip.text)
        tip._hide()
        self.assertIsNone(tip._tip)

    @unittest.skipUnless(HAS_SCIPY, "scipy not installed")
    def test_simulator_fit_recovers_size_and_amplitude(self):
        lam = np.arange(350.0, 801.0, 5.0)
        truth = 1.7 * lsm.vesicle_amplitude(lam, 130.0, 0.2, 1, 10.0 / 40, 1.0)      # 1.7x the lipid assumed
        win = self.open_simulator(lam, truth)
        v = win._sim_vars
        v['dil'].set('40')
        win._sim_fit(True)
        self.assertAlmostEqual(float(v['d'].get()), 130.0, delta=3.0)
        self.assertAlmostEqual(float(v['k'].get()), 1.7, delta=0.05)


if __name__ == "__main__":
    unittest.main()
