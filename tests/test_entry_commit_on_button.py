"""A value typed in an entry (without Enter) is used when a button is clicked right after, with real Tk events.

A Tk button does not take the keyboard focus on click, so the entry's <FocusOut> never fired before the command ran.

    python -m unittest discover -s tests -v
"""
import importlib.util
import os
import tempfile
import tkinter as tk
import unittest
from unittest import mock

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("lsm_focus", os.path.join(HERE, "..", "LabSpectrumManager.pyw"))
lsm = importlib.util.module_from_spec(spec)
spec.loader.exec_module(lsm)


class TestTrimEntryThenApply(unittest.TestCase):
    def setUp(self):
        try:
            self.root = lsm.TkinterDnD.Tk() if lsm.HAS_DND else tk.Tk()
        except tk.TclError as e:
            raise unittest.SkipTest("Tk not available: %s" % e)
        self.addCleanup(self.root.destroy)
        self.root.geometry("+3000+3000")                  # mapped (focus needs it) but off screen
        with mock.patch.object(lsm.messagebox, "showwarning"), mock.patch.object(lsm.messagebox, "askokcancel"):
            self.app = lsm.LabSpectrumManager(self.root)
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        path = os.path.join(tmp.name, "s.csv")
        with open(path, "w") as f:
            f.write("Wavelength [nm],S\n" + "".join(f"{w},{0.1 + w / 1000:.4f}\n" for w in range(250, 400)))
        self.app.processa_file(path)
        self.app.aggiorna_vista()
        self.root.update()

    def click(self, widget):
        x, y = widget.winfo_width() // 2, widget.winfo_height() // 2
        widget.event_generate('<Enter>', x=x, y=y)
        widget.event_generate('<ButtonPress-1>', x=x, y=y)
        self.root.update()
        widget.event_generate('<ButtonRelease-1>', x=x, y=y)
        self.root.update()

    def test_typed_trim_limit_is_used_without_moving_the_cursor(self):
        self.app.apri_trim()
        self.root.update()
        apply_btn = [w for w in self.app.f_trim.winfo_children()[-1].winfo_children() if w.cget('text') == 'Apply'][0]
        self.app._tr_entry_lo.focus_force()
        self.root.update()
        self.app._tr_entry_lo.delete(0, 'end')
        self.app._tr_entry_lo.insert(0, "300")            # no Enter, no click elsewhere
        self.click(apply_btn)
        self.assertEqual(self.app.spectra['S_trim']['df'].index.min(), 300)


if __name__ == "__main__":
    unittest.main()
