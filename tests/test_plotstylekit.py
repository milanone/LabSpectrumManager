"""How LabSpectrumManager works with and without PlotStyleKit (shared Origin style and figure editor).

The program carries the bootstrap of PlotStyleKit's integration module: when it cannot find PlotStyleKit it asks for
the folder, and without it every figure action has a fallback. The dialogs are replaced by answers; the application
runs on a hidden Tk window with a synthetic CSV spectrum.

    python -m unittest discover -s tests -v
"""
import importlib.util
import json
import os
import pickle
import tempfile
import tkinter as tk
import unittest
from unittest import mock

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("lsm_psk", os.path.join(HERE, "..", "LabSpectrumManager.pyw"))
lsm = importlib.util.module_from_spec(spec)
spec.loader.exec_module(lsm)

KIT = os.path.normpath(os.path.join(HERE, "..", "..", "PlotStyleKit"))


def write_csv(path):
    wl = np.arange(250, 400)
    with open(path, "w", newline="") as f:
        f.write("Wavelength [nm],Sample X1\n")
        for w in wl:
            f.write(f"{w},{0.3 + 0.2 * np.exp(-((w - 320) / 20.0) ** 2):.4f}\n")


@unittest.skipUnless(os.path.isfile(os.path.join(KIT, "plotstylekit.py")), "needs the PlotStyleKit folder next to this project")
class TestWithAndWithoutPlotStyleKit(unittest.TestCase):
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
        self.saved = (lsm.psk, lsm.PSK_PROBLEM, lsm.PSK_CONFIG)
        lsm.PSK_CONFIG = os.path.join(self.tmp.name, "psk.json")           # never touch the real ~/.plotstylekit.json
        self.addCleanup(self._restore)
        self.app = lsm.LabSpectrumManager(self.root)
        write_csv(self.path("s.csv"))
        self.app.processa_file(self.path("s.csv"))
        self.app.aggiorna_vista()

    def _restore(self):
        lsm.psk, lsm.PSK_PROBLEM, lsm.PSK_CONFIG = self.saved
        for w in self.root.winfo_children():
            if isinstance(w, tk.Toplevel):
                w.destroy()

    def path(self, name):
        return os.path.join(self.tmp.name, name)

    def test_the_kit_is_in_use(self):
        self.assertIsNotNone(lsm.psk)
        self.assertEqual(len(self.app.spectra), 1)

    def test_saved_figure_has_the_publication_size_and_no_cursor(self):
        n_lines = len(self.app.ax.lines)
        self.app.v_line = self.app.ax.axvline(300)                                # the interactive cursor
        with mock.patch.object(lsm.filedialog, 'asksaveasfilename', return_value=self.path('f.fig.pickle')), \
                mock.patch.object(lsm.messagebox, 'showinfo'):
            self.app.salva_figura_pickle()
        with open(self.path('f.fig.pickle'), 'rb') as f:
            fig = pickle.load(f)
        self.assertAlmostEqual(fig.get_size_inches()[0] * 2.54, 9.0, places=3)
        self.assertEqual(len(fig.axes[0].lines), n_lines)                         # what the plot has, not the cursor

    def test_without_the_kit_the_figure_is_saved_as_it_is(self):
        lsm.psk = None
        live = tuple(self.app.fig.get_size_inches())
        with mock.patch.object(lsm.filedialog, 'asksaveasfilename', return_value=self.path('plain.fig.pickle')), \
                mock.patch.object(lsm.messagebox, 'showinfo') as info, \
                mock.patch.object(lsm.messagebox, 'showerror') as err:
            self.app.salva_figura_pickle()
        self.assertFalse(err.called, err.call_args)
        self.assertIn('non ristilizzata', info.call_args[0][1])
        with open(self.path('plain.fig.pickle'), 'rb') as f:
            self.assertEqual(tuple(pickle.load(f).get_size_inches()), live)

    def test_without_the_kit_edit_figure_offers_to_locate_it(self):
        lsm.psk = None
        with mock.patch.object(lsm.messagebox, 'askyesno', return_value=False) as ask:
            self.app.apri_editor_figura()
        self.assertTrue(ask.called)
        self.assertIsNone(lsm.psk)
        self.assertFalse(any(isinstance(w, tk.Toplevel) for w in self.root.winfo_children()))

    def test_edit_figure_opens_the_editor_with_the_kit(self):
        with mock.patch.object(lsm.messagebox, 'showerror') as err:
            self.app.apri_editor_figura()
        self.assertFalse(err.called, err.call_args)
        self.assertTrue(any(isinstance(w, tk.Toplevel) for w in self.root.winfo_children()))

    def test_cancelling_the_warning_is_remembered(self):
        lsm.psk = None
        with mock.patch.object(lsm.messagebox, 'askokcancel', return_value=False) as ask:
            self.app._offer_plotstylekit("PlotStyleKit was not found")
        self.assertEqual(ask.call_count, 1)
        self.assertIsNone(lsm.psk)
        with open(lsm.PSK_CONFIG, encoding='utf-8') as f:
            self.assertTrue(json.load(f)['declined'])

    def test_choosing_the_folder_starts_the_kit_and_remembers_it(self):
        lsm.psk = None
        with mock.patch.object(lsm.messagebox, 'askokcancel', return_value=True), \
                mock.patch.object(lsm.filedialog, 'askdirectory', return_value=KIT), \
                mock.patch.object(lsm.messagebox, 'showinfo'):
            self.app._offer_plotstylekit("PlotStyleKit was not found")
        self.assertIsNotNone(lsm.psk)
        with open(lsm.PSK_CONFIG, encoding='utf-8') as f:
            cfg = json.load(f)
        self.assertEqual(os.path.normpath(cfg['path']), KIT)
        self.assertFalse(cfg['declined'])

    def test_a_wrong_folder_is_refused_with_the_reason(self):
        lsm.psk = None
        empty = self.path('empty')
        os.makedirs(empty)
        texts, answers = [], iter([True, False])
        with mock.patch.object(lsm.messagebox, 'askokcancel',
                               side_effect=lambda *a, **k: (texts.append(a[1]) or next(answers))), \
                mock.patch.object(lsm.filedialog, 'askdirectory', return_value=empty):
            self.app._offer_plotstylekit("PlotStyleKit was not found")
        self.assertEqual(len(texts), 2)
        self.assertIn("is not in", texts[1])
        self.assertIsNone(lsm.psk)


if __name__ == "__main__":
    unittest.main()
