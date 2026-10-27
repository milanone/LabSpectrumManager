"""The scattering correction reports its equation: the text, the panel, the clipboard and the metadata of the result.

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
spec = importlib.util.spec_from_file_location("lsm_scatter", os.path.join(HERE, "..", "LabSpectrumManager.pyw"))
lsm = importlib.util.module_from_spec(spec)
spec.loader.exec_module(lsm)


class TestEquationText(unittest.TestCase):
    def test_signs_are_folded_into_the_operators(self):
        t = lsm.scatter_equation(0.0123, 2.31, 1.2e-05, -0.021)
        self.assertEqual(t, "Scatter(λ) = 0.0123·(λ/1000)^(−2.31) − 1.2e-05·λ − 0.021   (λ in nm)")
        t = lsm.scatter_equation(0.5, 4.0, -3e-06, 0.1)
        self.assertEqual(t, "Scatter(λ) = 0.5·(λ/1000)^(−4) + 3e-06·λ + 0.1   (λ in nm)")

    def test_the_text_evaluates_to_the_subtracted_term(self):
        VS, P, M, OFF = 0.02, 2.5, 4e-06, -0.01
        x = np.array([300.0, 450.0, 700.0])
        rhs = lsm.scatter_equation(VS, P, M, OFF).split("=", 1)[1].split("(λ in nm)")[0]
        rhs = rhs.replace("·", "*").replace("^", "**").replace("−", "-")
        got = [eval(rhs.replace("λ", f"({v})")) for v in x]
        np.testing.assert_allclose(got, VS * (x / 1000) ** -P - M * x + OFF, rtol=1e-5)


class TestPanel(unittest.TestCase):
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
        wl = np.arange(250, 700)
        y = 0.05 + 0.4 * np.exp(-((wl - 320) / 25.0) ** 2) + 0.03 * (wl / 1000.0) ** -3
        path = os.path.join(self.tmp.name, "s.csv")
        with open(path, "w", newline="") as f:
            f.write("Wavelength [nm],Sample X1\n")
            for w, v in zip(wl, y):
                f.write(f"{w},{v:.5f}\n")
        self.app.processa_file(path)
        self.app.aggiorna_vista()

    def open_panel(self):
        name = next(iter(self.app.spectra))
        self.app.apri_correzione_scattering()          # aggiorna_vista() selected the only spectrum
        return name

    def test_equation_follows_the_parameters_and_ends_up_in_the_metadata(self):
        name = self.open_panel()
        self.app._linked_set(self.app._sc_VS, 0.0123)
        self.app._linked_set(self.app._sc_P, 2.5)
        shown = self.app._sc_eq_label.cget("text")
        self.assertIn("0.0123·(λ/1000)^(−2.5)", shown)
        self.app._sc_applica()
        info = self.app.spectra[name + "_corr"]["info"]
        self.assertIn("0.0123·(λ/1000)^(−2.5)", info["Scatter equation"])

    def test_auto_fit_reports_the_window_and_the_residual(self):
        self.open_panel()
        self.app._sc_auto()
        text = self.app._sc_fit_label.cget("text")
        self.assertIn("rms residual", text)
        self.assertIn("P optimised", text)
        self.assertIn("(λ/1000)^(−", self.app._sc_eq_label.cget("text"))

    def test_copy_puts_form_and_parameters_on_the_clipboard(self):
        self.open_panel()
        self.app._linked_set(self.app._sc_VS, 0.02)
        self.app._sc_copia_equazione()
        clip = self.root.clipboard_get()
        self.assertIn(lsm.SCATTER_FORM, clip)
        self.assertIn("VS = 0.02", clip)


if __name__ == "__main__":
    unittest.main()
