"""Help > Check for updates: version comparison, the GitHub queries (mocked, no network) and the dialogs."""
import importlib.util
import json
import os
import tempfile
import tkinter as tk
import unittest
import urllib.error
from unittest import mock

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("lsm", os.path.join(os.path.dirname(HERE), "LabSpectrumManager.pyw"))
lsm = importlib.util.module_from_spec(spec)
spec.loader.exec_module(lsm)


class TestVersion(unittest.TestCase):
    def test_parse_and_compare(self):
        self.assertEqual(lsm.parse_version("v1.2.3"), (1, 2, 3))
        self.assertEqual(lsm.parse_version("1.10"), (1, 10))
        self.assertEqual(lsm.parse_version("v2.0.0-rc1"), (2, 0, 0))
        self.assertGreater(lsm.parse_version("v1.10.0"), lsm.parse_version("v1.9.9"))     # numbers, not text
        self.assertEqual(lsm.parse_version("nonsense"), (0,))

    def test_the_program_has_a_version(self):
        self.assertTrue(lsm.parse_version(lsm.__version__) > (0,))


class TestFetch(unittest.TestCase):
    def test_latest_release(self):
        body = json.dumps({"tag_name": "v1.1.0", "body": "notes", "html_url": "https://example/r"}).encode()
        with mock.patch.object(lsm, "_github_get", return_value=body):
            rel = lsm.fetch_latest_release()
        self.assertEqual(rel, {"tag": "v1.1.0", "notes": "notes", "url": "https://example/r"})

    def test_no_release_falls_back_to_the_newest_tag(self):
        calls = []

        def fake(url, timeout=10):
            calls.append(url)
            if url.endswith("/releases/latest"):
                raise urllib.error.HTTPError(url, 404, "Not Found", {}, None)
            return json.dumps([{"name": "v1.9.0"}, {"name": "v1.10.0"}]).encode()
        with mock.patch.object(lsm, "_github_get", side_effect=fake):
            rel = lsm.fetch_latest_release()
        self.assertEqual(rel["tag"], "v1.10.0")                                          # not the textual maximum
        self.assertEqual(len(calls), 2)

    def test_nothing_published(self):
        def fake(url, timeout=10):
            if url.endswith("/releases/latest"):
                raise urllib.error.HTTPError(url, 404, "Not Found", {}, None)
            return b"[]"
        with mock.patch.object(lsm, "_github_get", side_effect=fake):
            self.assertIsNone(lsm.fetch_latest_release())

    def test_other_errors_propagate(self):
        err = urllib.error.HTTPError("u", 403, "rate limit", {}, None)
        with mock.patch.object(lsm, "_github_get", side_effect=err):
            with self.assertRaises(urllib.error.HTTPError):
                lsm.fetch_latest_release()


class TestDownload(unittest.TestCase):
    GOOD = b"class LabSpectrumManager:\n    pass\n"

    def test_saves_under_a_new_name_and_never_overwrites(self):
        with tempfile.TemporaryDirectory() as d, mock.patch.object(lsm, "_github_get", return_value=self.GOOD):
            p1 = lsm.download_release("v1.1.0", d)
            p2 = lsm.download_release("v1.1.0", d)
            self.assertEqual(os.path.basename(p1), "LabSpectrumManager_v1.1.0.pyw")
            self.assertNotEqual(p1, p2)
            with open(p1, "rb") as f:
                self.assertEqual(f.read(), self.GOOD)

    def test_a_file_that_does_not_compile_or_is_not_the_program_is_refused(self):
        for bad in (b"def broken(:\n", b"print('hello')\n"):
            with tempfile.TemporaryDirectory() as d, mock.patch.object(lsm, "_github_get", return_value=bad):
                with self.assertRaises((SyntaxError, ValueError)):
                    lsm.download_release("v1.1.0", d)
                self.assertEqual(os.listdir(d), [])                                       # nothing was written


class TestDialogs(unittest.TestCase):
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
        with mock.patch.object(lsm.messagebox, "showwarning"):
            self.app = lsm.LabSpectrumManager(self.root)

    def test_the_menu_has_the_entry(self):
        self.assertEqual(self.app.help_menu.entrycget(0, "label"), "Check for updates...")

    def test_up_to_date_error_and_nothing_published_are_messages(self):
        with mock.patch.object(lsm.messagebox, "showinfo") as info, mock.patch.object(lsm.messagebox, "showwarning") as warn:
            self.app._update_result({"tag": "v" + lsm.__version__, "notes": "", "url": "u"}, None)
            self.assertIn("up to date", info.call_args[0][1])
            self.app._update_result(None, None)
            self.assertIn("No published version", info.call_args[0][1])
            self.app._update_result(None, OSError("no network"))
            self.assertIn("no network", warn.call_args[0][1])

    def test_a_newer_version_opens_the_dialog_and_download_asks_for_a_folder(self):
        before = set(self.root.winfo_children())
        self.app._update_result({"tag": "v99.0.0", "notes": "new things", "url": "https://example/r"}, None)
        win = [w for w in self.root.winfo_children() if w not in before and isinstance(w, tk.Toplevel)][0]
        self.assertEqual(win._update_rel["tag"], "v99.0.0")
        with tempfile.TemporaryDirectory() as d, mock.patch.object(lsm.filedialog, "askdirectory", return_value=d), \
                mock.patch.object(lsm, "_github_get", return_value=TestDownload.GOOD), \
                mock.patch.object(lsm.messagebox, "showinfo") as info:
            self.app._update_download("v99.0.0", win)
            self.assertTrue(os.path.exists(os.path.join(d, "LabSpectrumManager_v99.0.0.pyw")))
            self.assertIn("running program was not changed", info.call_args[0][1])
        win.destroy()


if __name__ == "__main__":
    unittest.main()
