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

    def test_the_exponent_has_no_upper_limit(self):
        y = 0.5 * (self.x / 1000.0) ** -7.5
        vs, p, m, rms = lsm.fit_power_law(self.x, y)
        self.assertAlmostEqual(p, 7.5, places=2)
        self.assertAlmostEqual(vs, 0.5, places=3)
        self.assertLess(rms, 1e-4)

    def test_under_keeps_the_trace_below_the_points_and_gives_way_with_vs(self):
        y = 0.5 * (self.x / 1000.0) ** -3.0
        y[40] -= 0.03                                                     # one point below the power law
        free = lsm.fit_power_law(self.x, y)
        vs, p, m, rms = lsm.fit_power_law(self.x, y, under=(self.x, y))
        trace = vs * (self.x / 1000.0) ** -p - m * self.x
        self.assertLessEqual(float(np.max(trace - y)), 1e-12)
        self.assertLess(vs, free[0])                                      # VS gave way
        self.assertLess(float(np.min(y - trace)), 1e-3)                  # and the trace touches the data

    def test_under_leaves_an_exact_power_law_alone(self):
        y = 0.5 * (self.x / 1000.0) ** -3.0
        vs, p, m, rms = lsm.fit_power_law(self.x, y, under=(self.x, y))
        self.assertAlmostEqual(p, 3.0, places=3)
        self.assertAlmostEqual(vs, 0.5, places=3)

    def test_under_holds_after_rounding_also_with_a_free_slope_and_a_fixed_exponent(self):
        rng = np.random.default_rng(3)
        y = 0.05 * (self.x / 1000.0) ** -4.0 + 0.1 * np.exp(-0.5 * ((self.x - 520) / 25) ** 2) + rng.normal(0, 0.001, self.x.size)
        for kw in ({'m': None}, {'P': 3.3}, {'m': 0.0}):
            vs, p, m, rms = lsm.fit_power_law(self.x, y, under=(self.x, y), **kw)
            trace = vs * (self.x / 1000.0) ** -p - m * self.x
            self.assertLessEqual(float(np.max(trace - y)), 1e-12, msg=str(kw))

    def test_under_also_applies_to_points_outside_the_fitted_window(self):
        y = 0.5 * (self.x / 1000.0) ** -3.0
        y[-1] -= 0.05                                                     # the last point, left out of the window
        win = self.x < 780
        vs, p, m, rms = lsm.fit_power_law(self.x[win], y[win], under=(self.x, y))
        trace = vs * (self.x / 1000.0) ** -p - m * self.x
        self.assertLessEqual(float(np.max(trace - y)), 1e-12)

    def test_under_leaves_out_the_points_at_or_below_zero(self):
        y = 0.5 * (self.x / 1000.0) ** -3.0
        y[-30:] = -0.001                                                  # a blank-corrected tail below zero
        for kw in ({'m': 0.0}, {'m': None}):
            vs, p, m, rms = lsm.fit_power_law(self.x, y, under=(self.x, y), **kw)
            self.assertGreater(vs, 0.0, msg=str(kw))                      # no collapse to VS = 0
            pos = y > 0
            trace = vs * (self.x / 1000.0) ** -p - m * self.x
            self.assertLessEqual(float(np.max((trace - y)[pos])), 1e-12, msg=str(kw))

    def test_a_limit_can_still_be_asked_for(self):
        y = 0.5 * (self.x / 1000.0) ** -7.5
        self.assertLessEqual(lsm.fit_power_law(self.x, y, p_hi=5.0)[1], 5.0)

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
class TestSolidSphere(unittest.TestCase):
    lam = np.array([350.0, 400.0, 500.0, 600.0, 800.0])

    def test_reference_values_of_a_100_nm_latex(self):
        """0.005% w/v, 100 nm, monodisperse, 1 cm (the values printed on the dilution sheet, n(589) = 1.596)."""
        a = lsm.solid_sphere_amplitude(self.lam, 100.0, 0.0, 0.005, n589=1.596)
        np.testing.assert_allclose(a, [0.2246, 0.1356, 0.0566, 0.0274, 0.0086], rtol=0.01)

    def test_proportional_to_concentration_and_path(self):
        a = lsm.solid_sphere_amplitude(self.lam, 100.0, 0.05, 0.01)
        np.testing.assert_allclose(lsm.solid_sphere_amplitude(self.lam, 100.0, 0.05, 0.04, path_cm=0.5), 2.0 * a,
                                   rtol=1e-9)

    def test_matches_the_textbook_mie_sum_for_a_homogeneous_sphere(self):
        from scipy.special import spherical_jn, spherical_yn
        lam = np.array([500.0])
        d, n589, rho, pct = 100.0, 1.592, 1.055, 0.005
        nm = lsm._water_index(lam)[0]
        n_p = n589 + (lsm._ps_dispersion(500.0) - lsm._ps_dispersion(589.0))
        m_rel = n_p / nm
        x = 2 * np.pi * nm / 500.0 * d / 2
        s = 0.0
        for n in range(1, int(x + 4 * x ** (1 / 3) + 2) + 1):
            def ric(z):
                j, y = spherical_jn(n, z), spherical_yn(n, z)
                dj, dy = spherical_jn(n, z, derivative=True), spherical_yn(n, z, derivative=True)
                return z * j, j + z * dj, z * (j - 1j * y), (j - 1j * y) + z * (dj - 1j * dy)
            pm, dm, _, _ = ric(m_rel * x)
            p, dp, xi, dxi = ric(x)
            a = (m_rel * pm * dp - p * dm) / (m_rel * pm * dxi - xi * dm)
            b = (pm * dp - m_rel * p * dm) / (pm * dxi - m_rel * xi * dm)
            s += (2 * n + 1) * (abs(a) ** 2 + abs(b) ** 2)
        c_sca = 2 / x ** 2 * s * np.pi * (d / 2) ** 2                                    # nm^2
        volume = np.pi / 6 * d ** 3 * 1e-21
        ref = (pct / 100 / rho / volume) * c_sca * 1e-14 / np.log(10)
        got = lsm.solid_sphere_amplitude(lam, d, 0.0, pct, 1.0, n589, rho, acc_deg=0.0)[0]
        self.assertAlmostEqual(got / ref, 1.0, delta=0.01)

    def test_the_apparent_P_of_a_100_nm_polystyrene_sphere_is_near_four(self):
        lam = np.arange(350.0, 801.0, 10.0)
        P = lsm.fit_power_law(lam, lsm.solid_sphere_amplitude(lam, 100.0, 0.0, 0.01), m=0.0)[1]
        self.assertTrue(3.7 < P < 4.1, P)


@unittest.skipUnless(HAS_SCIPY, "scipy not installed")
class TestAmplitude(unittest.TestCase):
    lam = np.array([500.0])

    def test_unilamellar_reference_value(self):
        """Value computed while developing the model (100 nm, 10 mg/mL POPC, 0.1 cm cuvette, shell index excess 0.155;
        it was 0.0606 with the earlier default 0.12, and the amplitude goes as the excess squared)."""
        a = lsm.vesicle_amplitude(self.lam, 100.0, 0.2, 1, 10.0, 0.1)[0]
        self.assertAlmostEqual(a, 0.1007, delta=0.002)
        a_old = lsm.vesicle_amplitude(self.lam, 100.0, 0.2, 1, 10.0, 0.1, dn=0.12)[0]
        self.assertAlmostEqual(a / a_old, (0.155 / 0.12) ** 2, delta=0.12)

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
        self.assertEqual((mx['mw'], mx['area'], mx['t'], mx['dn'], mx['peg_film']), (760.08, 0.65, 4.5, 0.155, 0.0))

    def test_85_10_5_mixture(self):
        mx = lsm.lipid_mixture(10, 5)
        self.assertAlmostEqual(mx['mw'], 825.0, delta=0.5)            # 0.85*760.08 + 0.10*386.65 + 0.05*2805.5
        self.assertAlmostEqual(mx['area'], 0.5915, places=3)
        self.assertAlmostEqual(mx['t'], 4.9, places=6)
        self.assertAlmostEqual(mx['peg_film'], 0.234, delta=0.01)     # nm of pure PEG per leaflet
        self.assertGreater(mx['dn'], 0.155 + 0.01)                     # the PEG film adds polarisability

    def test_fractions_above_100_percent_are_refused(self):
        with self.assertRaises(ValueError):
            lsm.lipid_mixture(60, 50)

    def test_no_peg_leaflets_removes_the_film_term(self):
        a = lsm.lipid_mixture(10, 5, peg_leaflets=0)
        b = lsm.lipid_mixture(10, 5, peg_leaflets=2)
        self.assertAlmostEqual(a['dn'], 0.155, places=9)
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

    def test_auto_follows_an_exponent_above_5_and_the_slider_grows(self):
        x = self.app._sc_x
        self.app._sc_y = 0.5 * (x / 1000.0) ** -7.0
        self.app._sc_OFF_lock.set(True)
        self.app._sc_auto()
        self.assertAlmostEqual(self.app._sc_P['value'].get(), 7.0, places=2)
        self.assertGreaterEqual(self.app._sc_P['hi'], 7.0)             # the slider range was extended

    def test_the_general_auto_fit_is_not_limited_to_5_either(self):
        x = self.app._sc_x
        self.app._sc_y = 0.0015 * (x / 1000.0) ** -6.0 + 0.002        # a realistic amplitude: 0.8 at 350 nm
        self.app._sc_auto()
        self.assertAlmostEqual(self.app._sc_P['value'].get(), 6.0, places=1)
        self.assertGreaterEqual(self.app._sc_P['hi'], 6.0)

    def test_vs_follows_a_very_large_exponent_in_both_fits(self):
        x = self.app._sc_x
        for lock, P_true, off in ((True, 20, 0.0), (False, 12, 0.0), (False, 20, 0.003)):
            vs = 0.35 ** P_true                                           # an absorbance of 1 at 350 nm
            self.app._sc_y = vs * (x / 1000.0) ** -P_true + off
            self.app._linked_set(self.app._sc_P, 2.0)
            self.app._sc_OFF_lock.set(lock)
            self.app._sc_auto()
            self.assertAlmostEqual(self.app._sc_P['value'].get(), P_true, delta=0.15 if lock and off else 0.02)
            self.assertAlmostEqual(self.app._sc_VS['value'].get() / vs, 1.0, delta=0.3 if off else 0.01)
            # the slider range follows VS, and the box holds it with enough digits (a re-read must not zero it)
            self.assertLess(self.app._sc_VS['hi'], 1e-3)
            self.assertAlmostEqual(float(self.app._sc_VS['entry'].get()) / self.app._sc_VS['value'].get(), 1.0, places=4)

    def test_reopening_the_panel_resets_the_vs_range(self):
        x = self.app._sc_x
        self.app._sc_y = 0.35 ** 12 * (x / 1000.0) ** -12.0
        self.app._sc_OFF_lock.set(True)
        self.app._sc_auto()
        self.assertLess(self.app._sc_VS['hi'], 1e-3)
        self.app._sc_annulla()
        self.app.apri_correzione_scattering()
        self.assertEqual(self.app._sc_VS['hi'], 1.0)

    def test_the_general_auto_fit_keeps_the_trace_under_the_spectrum_after_rounding(self):
        x = self.app._sc_x
        rng = np.random.default_rng(1)
        cases = (0.07 * (x / 1000.0) ** -2.6 + 0.01 + rng.normal(0, 0.002, x.size),
                 0.05 * (x / 1000.0) ** -4 + 0.15 * np.exp(-0.5 * ((x - 520) / 25) ** 2) + 0.01 + rng.normal(0, 0.001, x.size),
                 0.35 ** 10 * (x / 1000.0) ** -10 + 0.004 + rng.normal(0, 0.0005, x.size))
        for y in cases:
            self.app._sc_y = y
            self.app._linked_set(self.app._sc_P, 2.0)
            self.app._sc_OFF_lock.set(False)
            self.app._sc_auto()
            VS, P, M, OFF = (self.app._sc_VS['value'].get(), self.app._sc_P['value'].get(),
                             self.app._sc_M['value'].get(), self.app._sc_OFF['value'].get())
            trace = VS * (x / 1000.0) ** -P - M * x + OFF                  # with the values as shown in the panel
            self.assertLessEqual(float(np.max(trace - y)), 1e-9)

    def test_keep_under_data_flag_is_only_available_with_the_offset_fixed_and_forces_the_trace_under(self):
        x = self.app._sc_x
        self.assertEqual(str(self.app._sc_under_cb.cget("state")), "disabled")
        self.app._sc_OFF_lock.set(True)
        self.assertEqual(str(self.app._sc_under_cb.cget("state")), "normal")
        rng = np.random.default_rng(2)
        y = 0.07 * (x / 1000.0) ** -2.6 + rng.normal(0, 0.002, x.size)
        self.app._sc_y = y
        self.app._sc_auto()
        VS, P = self.app._sc_VS['value'].get(), self.app._sc_P['value'].get()
        self.assertGreater(float(np.max(VS * (x / 1000.0) ** -P - y)), 0.0)  # plain least squares crosses the noise
        self.app._sc_UNDER.set(True)
        self.app._sc_auto()
        VS, P, M = self.app._sc_VS['value'].get(), self.app._sc_P['value'].get(), self.app._sc_M['value'].get()
        self.assertLessEqual(float(np.max(VS * (x / 1000.0) ** -P - M * x - y)), 1e-9)      # the panel's own values
        self.assertEqual(self.app._sc_OFF['value'].get(), 0.0)
        self.assertIn("kept under the data", self.app._sc_fit_label.cget("text"))

    def test_keep_under_data_starts_also_when_some_points_are_not_positive(self):
        x = self.app._sc_x
        y = 0.07 * (x / 1000.0) ** -2.6
        y[-40:] = -0.0005                                                  # a negative tail, as after a blank subtraction
        self.app._sc_y = y
        self.app._sc_OFF_lock.set(True)
        self.app._sc_UNDER.set(True)
        self.app._sc_auto()
        VS, P = self.app._sc_VS['value'].get(), self.app._sc_P['value'].get()
        self.assertGreater(VS, 0.0)
        self.assertLessEqual(float(np.max((VS * (x / 1000.0) ** -P - y)[y > 0])), 1e-9)
        self.assertIn("40 points at or below zero", self.app._sc_fit_label.cget("text"))

    def test_ticking_the_box_zeroes_the_offset(self):
        self.app._linked_set(self.app._sc_OFF, 0.05)
        self.app._sc_OFF_lock.set(True)
        self.app._sc_OFF_lock.get() and self.app._linked_set(self.app._sc_OFF, 0.0)   # what the checkbox command does
        self.assertEqual(self.app._sc_OFF['value'].get(), 0.0)

class TestSimulator(unittest.TestCase):
    """The spectrum simulator is a tool of its own: it works on the selected spectrum, with no scattering panel open."""

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

    def load(self, x, y):
        """Load (x, y) as the only spectrum of the app, so that it is the selected one."""
        path = os.path.join(self.tmp.name, "sim.csv")
        with open(path, "w", newline="") as f:
            f.write("Wavelength [nm],Sample X1\n")
            for w, v in zip(x, y):
                f.write(f"{w},{v:.9f}\n")
        self.app.processa_file(path)
        self.app.aggiorna_vista()

    def open_simulator(self, x, y):
        """Open the simulator on the spectrum (x, y); returns the window."""
        self.load(x, y)
        before = set(self.root.winfo_children())
        self.app.open_simulator()
        win = [w for w in self.root.winfo_children() if w not in before and isinstance(w, tk.Toplevel)][0]
        self.addCleanup(win.destroy)
        return win

    def test_it_is_independent_of_the_scattering_panel(self):
        lam = np.arange(350.0, 801.0, 5.0)
        win = self.open_simulator(lam, 0.07 * (lam / 1000.0) ** -2.6)
        self.assertIsNone(self.app._sc_name)                           # no scattering panel was opened
        self.assertEqual((win._sim_vars['lo'].get(), win._sim_vars['hi'].get()), ('350', '800'))

    def test_no_selection_or_ftir_is_refused(self):
        with mock.patch.object(lsm.messagebox, "showwarning") as warn:
            before = set(self.root.winfo_children())
            self.app.open_simulator()                                  # nothing loaded
        warn.assert_called_once()
        self.assertEqual(set(self.root.winfo_children()), before)

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
        # label + entry per field, the sample-type label and combobox, 3 buttons
        self.assertEqual(len(win._sim_tips), 2 * len(win._sim_vars) + 2 + 3)
        tip = next(t for t in win._sim_tips if t.text == lsm.SIM_TIPS['sg'])
        tip._show()
        self.assertIsNotNone(tip._tip)
        self.assertIn("PDI", tip.text)
        tip._hide()
        self.assertIsNone(tip._tip)

    @unittest.skipUnless(HAS_SCIPY, "scipy not installed")
    def test_simulator_solid_sphere_mode_overlays_the_latex_theory(self):
        lam = np.arange(350.0, 801.0, 10.0)
        truth = lsm.solid_sphere_amplitude(lam, 100.0, 0.05, 0.005)             # 0.005% w/v, 1 cm, the default n and density
        win = self.open_simulator(lam, truth)
        v = win._sim_vars
        win._sim_mode.set(lsm.SIM_MODES[1])
        self.assertEqual((v['c'].get(), v['sg'].get()), ('4', '0.05'))            # usual starting values of a latex
        v['dil'].set('800'); v['d'].set('100')                                    # 4% w/v / 800 = 0.005% w/v
        win._sim_update()
        text = [w.cget("text") for w in win.winfo_children() if isinstance(w, tk.Label)][-1]
        self.assertIn("Solid spheres in the cuvette 0.005% w/v", text)
        rms = float(text.split("rms residual")[1].split()[0])
        self.assertLess(rms, 1e-3 * truth.max())
        win._sim_mode.set(lsm.SIM_MODES[0])                                      # and back to the vesicles
        self.assertEqual(v['c'].get(), '10')

    @unittest.skipUnless(HAS_SCIPY, "scipy not installed")
    def test_simulator_fit_finds_the_diameter_of_a_latex(self):
        lam = np.arange(350.0, 801.0, 10.0)
        truth = lsm.solid_sphere_amplitude(lam, 125.0, 0.05, 0.005)
        win = self.open_simulator(lam, truth)
        v = win._sim_vars
        win._sim_mode.set(lsm.SIM_MODES[1])
        v['dil'].set('800')
        win._sim_fit(True)
        self.assertAlmostEqual(float(v['d'].get()), 125.0, delta=3.0)
        self.assertAlmostEqual(float(v['k'].get()), 1.0, delta=0.05)

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
