"""Test del lettore di export CSV BMG (micropiastre) di LabSpectrumManager.

I file sono sintetici, costruiti con la stessa struttura dell'export reale (vedi la sezione
"BMG plate reader" del README); non servono dati veri né una finestra Tk.

    python -m unittest discover -s tests -v
"""
import importlib.util
import os
import tempfile
import unittest
from unittest import mock

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location(
    "lsm", os.path.join(HERE, "..", "LabSpectrumManager.pyw"))
lsm = importlib.util.module_from_spec(spec)
spec.loader.exec_module(lsm)

WL = np.arange(350, 360, 1.0)          # 10 lunghezze d'onda
BLANKS = {"A01": 0.30, "B01": 0.40, "C01": 0.50}   # valori costanti per riga
SAMPLES = {"A02": 0.60, "B02": 0.70}


def fmt(v, dec="."):
    return f"{v:.3f}".replace(".", dec)


def bmg_text(corrected=False, dec=".", with_path=True):
    """Testo di un export BMG: grezzo, oppure con bianco già sottratto (righe Blank vuote)."""
    media = np.mean(list(BLANKS.values()))
    kind = ("Blank corrected based on Raw Data   (Abs Spectrum)#" if corrected
            else "Raw Data  (Abs Spectrum)#")
    path = r"Path: C:\Program Files (x86)\BMG\Omega\User\Data;" if with_path else ""
    out = [
        f"User: USER;{path}Test run no.: 7",
        "Test name: demo_350_360;Date: 03/07/2025;Time: 17:02:38",
        "",
        "ID1: demo experiment;",
        "Absorbance spectrum",
        "",
        kind,
        "",
        " Well; Content; " + "; ".join(["1 - 1"] * len(WL)) + "; ",
        ";Wavelength [nm];" + ";".join(f"{w:g}" for w in WL) + ";",
    ]
    for well, v in BLANKS.items():
        cells = [""] * len(WL) if corrected else [fmt(v, dec)] * len(WL)
        out.append(f"{well};Blank B;" + ";".join(cells) + ";")
    for well, v in SAMPLES.items():
        y = v - media if corrected else v
        out.append(f"{well};Sample {well};" + ";".join([fmt(y, dec)] * len(WL)) + ";")
    return "\r\n".join(out) + "\r\n"


class BmgReader(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.app = lsm.LabSpectrumManager.__new__(lsm.LabSpectrumManager)
        self.app.spectra = {}

    def write(self, name, text):
        p = os.path.join(self.tmp.name, name)
        with open(p, "w", encoding="latin-1", newline="") as f:
            f.write(text)
        return p

    def values(self, label):
        return self.app.spectra[label]["df"][label].to_numpy()

    def test_raw_with_blank_subtraction(self):
        p = self.write("run.CSV", bmg_text())
        self.app._leggi_csv_bmg(p, sottrai_bianco=True)
        self.assertEqual(sorted(self.app.spectra), ["run_A02", "run_B02"])
        media = np.mean(list(BLANKS.values()))
        np.testing.assert_allclose(self.values("run_A02"), SAMPLES["A02"] - media, atol=1e-3)
        info = self.app.spectra["run_A02"]["info"]
        self.assertEqual(info["Type"], "UV-Vis")
        self.assertEqual(info["Content"], "Sample A02")
        self.assertEqual(info["Experiment"], "demo experiment")
        self.assertEqual(info["Date"], "03/07/2025 17:02:38")
        self.assertIn("Raw Data", info["Data"])
        self.assertIn("Blank", info)

    def test_raw_without_subtraction_loads_blanks_too(self):
        p = self.write("run.CSV", bmg_text())
        self.app._leggi_csv_bmg(p, sottrai_bianco=False)
        self.assertEqual(len(self.app.spectra), 5)
        np.testing.assert_allclose(self.values("run_A01"), BLANKS["A01"], atol=1e-3)
        np.testing.assert_allclose(self.values("run_A02"), SAMPLES["A02"], atol=1e-3)

    def test_vendor_corrected_equals_own_subtraction(self):
        raw = self.write("raw.CSV", bmg_text())
        cor = self.write("cor.CSV", bmg_text(corrected=True))
        self.app._leggi_csv_bmg(raw, sottrai_bianco=True)
        self.app._leggi_csv_bmg(cor)      # nessuna domanda: i bianchi sono righe vuote
        self.assertEqual(len(self.app.spectra), 4)
        np.testing.assert_allclose(self.values("raw_B02"), self.values("cor_B02"), atol=1e-3)

    def test_decimal_comma(self):
        p = self.write("it.CSV", bmg_text(dec=","))
        self.app._leggi_csv_bmg(p, sottrai_bianco=False)
        np.testing.assert_allclose(self.values("it_B02"), SAMPLES["B02"], atol=1e-3)

    def test_wavelength_axis(self):
        p = self.write("run.CSV", bmg_text())
        self.app._leggi_csv_bmg(p, sottrai_bianco=False)
        np.testing.assert_allclose(self.app.spectra["run_A02"]["df"].index.to_numpy(), WL)

    def test_loading_twice_does_not_overwrite(self):
        p = self.write("run.CSV", bmg_text())
        self.app._leggi_csv_bmg(p, sottrai_bianco=False)
        self.app._leggi_csv_bmg(p, sottrai_bianco=False)
        self.assertEqual(len(self.app.spectra), 10)
        self.assertIn("run_A02_2", self.app.spectra)

    def test_missing_wavelength_row_is_an_error(self):
        p = self.write("bad.CSV", "User: X;Test run no.: 1\r\nA01;Sample;0.1;0.2\r\n")
        with self.assertRaises(ValueError):
            self.app._leggi_csv_bmg(p)

    def test_leggi_csv_dispatches_and_asks_about_blanks(self):
        p = self.write("run.CSV", bmg_text())
        with mock.patch.object(lsm.messagebox, "askyesno", return_value=True) as ask:
            self.app.leggi_csv(p)
        ask.assert_called_once()
        self.assertEqual(sorted(self.app.spectra), ["run_A02", "run_B02"])

    def test_plain_csv_is_not_taken_for_bmg(self):
        p = self.write("plain.csv", "nm,A\n400,0.1\n401,0.2\n402,0.3\n")
        self.app._UNIT_RE = lsm.LabSpectrumManager._UNIT_RE
        self.app.leggi_csv(p)
        # il CSV generico usa il nome del file come etichetta (la colonna "A" è solo un'unità)
        self.assertEqual(list(self.app.spectra), ["plain"])
        self.assertNotIn("Well", self.app.spectra["plain"]["info"])


if __name__ == "__main__":
    unittest.main()
