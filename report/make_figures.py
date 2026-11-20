"""Builds the figures and the numbers of report/size_from_turbidity_report.md from the example data
(`spettri di esempio/`, not in the repository because it holds real samples).

Run:  py report/make_figures.py          (writes report/figures/*.png and prints the numbers used in the report)
"""
import glob
import importlib.util
import os
import re

import matplotlib
import matplotlib.ticker
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.optimize import minimize_scalar

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = os.path.dirname(HERE)
DATA = os.path.join(BASE, "spettri di esempio")
OUT = os.path.join(HERE, "figures")
os.makedirs(OUT, exist_ok=True)


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


m = _load("lsm", os.path.join(BASE, "LabSpectrumManager.pyw"))
for cand in (os.path.join(BASE, "origin_style.py"), os.path.join(os.path.dirname(BASE), "PlotStyleKit", "origin_style.py")):
    if os.path.exists(cand):                      # same lookup order as the program itself
        _load("origin_style", cand).applica_rcparams()
        break
plt.rcParams.update({"font.size": 16, "axes.titlesize": 16, "axes.labelsize": 16, "xtick.labelsize": 16,
                     "ytick.labelsize": 16, "legend.fontsize": 14, "lines.linewidth": 2.0})
COL = plt.rcParams["axes.prop_cycle"].by_key()["color"]

# ---- lot data (certificate of analysis 1229748, 17 Sep 2012) -------------------------------------------------------
TEM_D, TEM_CV, SOLIDS, N589 = 130.0, 0.104, 4.2, 1.591
SIGMA = TEM_CV                                   # log-normal width ~ CV for a narrow population


def read_dsp(folder, pattern="*.dsp"):
    app = m.LabSpectrumManager.__new__(m.LabSpectrumManager)
    app.spectra = {}
    for f in sorted(glob.glob(os.path.join(DATA, folder, pattern))):
        app.leggi_dsp(f)
    return {n: (d["df"].index.to_numpy(float), d["df"].iloc[:, 0].to_numpy(float)) for n, d in app.spectra.items()}


def csv_count_rate(folder, name, num):
    c = glob.glob(os.path.join(DATA, folder, glob.escape(name) + "_%04d.csv" % (num + 2)))
    if c:
        for ln in open(c[0], encoding="latin-1"):
            if ln.startswith("Count Rate"):
                return float(ln.split(",")[1])
    return np.nan


def read_dls(folder):
    out = {}
    for f in sorted(glob.glob(os.path.join(DATA, folder, "*.nsz"))):
        nm = os.path.basename(f)[:-4]
        try:
            r = m.read_nsz(open(f, "rb").read())
        except ValueError:                       # the blank has no Z-average
            continue
        key, num = re.match(r"(.*)_(\d{4})$", nm).groups()
        r["rec"] = {"tau": r["tau"], "acf": r["acf"], "angle": r["angle"], "n": r["solvent_index"],
                    "visc": r["viscosity"], "T": r["temperature"]}
        r["cr"] = csv_count_rate(folder, key, int(num))
        out.setdefault(key, []).append(r)
    return out


def shape_k(lam, y, a):
    """Least-squares scale k of the model a onto y and the rms left (relative to the mean of y)."""
    k = float(y @ a / (a @ a))
    return k, float(np.sqrt(np.mean((y - k * a) ** 2)) / y.mean())


def mark(ax, letter):
    ax.text(-0.11, 1.02, letter, transform=ax.transAxes, fontsize=16, fontweight="bold", va="bottom")


lam = np.arange(400.0, 801.0, 10.0)
res = {}

# =====================================================================================================================
# 1. Latex UV against the Mie prediction from the certificate
# =====================================================================================================================
uv_lat = {k.split(" lattice")[0]: v for k, v in read_dsp("20261119 lattice 100 nm").items()}
dil = {"L1": 800, "L2": 400, "L3": 200}
lat = {}
for k, dl in dil.items():
    c = SOLIDS / dl
    Y = np.interp(lam, *uv_lat[k])
    A = m.solid_sphere_amplitude(lam, TEM_D, SIGMA, c, 1.0, n589=N589)
    kk, rms = shape_k(lam, Y, A)
    best = min(np.arange(100, 200, 1.0), key=lambda d: abs(float(Y.mean() / m.solid_sphere_amplitude(
        lam, d, SIGMA, c, 1.0, n589=N589).mean()) - 1))
    sel = lam >= 400
    P = m.fit_power_law(lam[sel], Y[sel], m=0.0)[1]
    lat[k] = dict(dil=dl, c=c, Y=Y, A=A, k=float(Y.mean() / A.mean()), rms=rms, d_k1=float(best), P=P,
                  A532=float(np.interp(532, *uv_lat[k])))
res["latex_uv"] = {k: {x: v for x, v in d.items() if x not in ("Y", "A")} for k, d in lat.items()}

fig, ax = plt.subplots(2, 1, figsize=(9, 11.5))
for i, k in enumerate(dil):
    ax[0].plot(lam, lat[k]["Y"], "-", color=COL[i], lw=1.6, label=f"{k} (1:{dil[k]}) measured")
    ax[0].plot(lam, lat[k]["A"], "--", color=COL[i], lw=1.0)
    ax[1].plot(lam, lat[k]["Y"] / lat[k]["A"], color=COL[i], lw=1.4, label=f"{k}: k = {lat[k]['k']:.2f}")
ax[0].plot([], [], "k--", lw=1.0, label="Mie, 130 nm, CV 10.4 %,\n4.2 % w/v, no free parameter")
ax[0].set_xlabel("Wavelength (nm)")
ax[0].set_ylabel("Absorbance (1 cm)")
ax[0].legend(fontsize=14)
ax[1].axhline(1, color="k", lw=0.8)
ax[1].set_ylim(0.8, 1.5)
ax[1].set_xlabel("Wavelength (nm)")
ax[1].set_ylabel("Measured / predicted")
ax[1].legend(fontsize=14, loc="lower right")
mark(ax[0], "a")
mark(ax[1], "b")
fig.tight_layout(pad=1.6, h_pad=2.4)
fig.savefig(os.path.join(OUT, "fig1_latex_uv_vs_mie.png"), dpi=150)
plt.close(fig)

# =====================================================================================================================
# 2. Latex DLS: Z-average against dilution and angle, size distributions, consistency of the dilutions
# =====================================================================================================================
D19 = "20261119 DLS lattice e liposomi"
dls19 = read_dls(D19)
LAT_KEYS = {800: "Polistirene 1a800 in NaCl 10 mM", 400: "Polistirene 1a400 in NaCl 10 mM",
            200: "Polistirene 1a200 in NaCl 10 mM", 100: "Polistirene 1a100 in NaCl 10 mM"}
ND = {800: 1.0, 400: 0.10, 200: 0.10}
I0 = (31692 + 31811) / 2.0
dls_lat = {}
for dl, key in LAT_KEYS.items():
    rs = dls19[key]
    dls_lat[dl] = dict(z=float(np.mean([r["z_average"] for r in rs])), pi=float(np.mean([r["pi"] for r in rs])),
                       cr=float(np.mean([r["cr"] for r in rs])), n=len(rs), angle=rs[0]["angle"],
                       t_before=float(np.mean([r["transmission_before"] for r in rs])), rs=rs)
rs90 = dls19["polistirene 1a100 in NaCl 10 mM 90 gradi"]
res["latex_dls"] = {dl: {x: v for x, v in d.items() if x != "rs"} for dl, d in dls_lat.items()}
res["latex_dls_1to100_90deg"] = dict(z=float(np.mean([r["z_average"] for r in rs90])),
                                      pi=float(np.mean([r["pi"] for r in rs90])), n=len(rs90))

# scaling between the two dilutions (expected x2 for the concentration)
tau = {dl: float(np.log(I0 / dls_lat[dl]["t_before"])) for dl in (800, 400, 200)}
cr_nd = {dl: dls_lat[dl]["cr"] / ND[dl] for dl in (800, 400, 200)}
ratios = {"UV absorbance at 532 nm": lat["L2"]["A532"] / lat["L1"]["A532"],
          "DLS extinction ln(I0/T)": tau[400] / tau[800],
          "DLS count rate / ND": cr_nd[400] / cr_nd[800]}
res["dilution_ratios_L2_over_L1"] = ratios
# count rate / extinction with the self-attenuation of the 1 x 0.4 cm cell (0.5 cm in, 0.2 cm out)
R = {dl: cr_nd[dl] * np.exp(tau[dl] * 0.7) / tau[dl] for dl in (800, 400, 200)}
res["R_count_rate_over_extinction"] = R

fig, axd = plt.subplot_mosaic([["a"], ["b"], ["c"], ["d"]], figsize=(9, 21))
ax = [axd["a"], axd["b"], axd["c"]]
ax[0].axhspan(TEM_D - 14, TEM_D + 14, color="0.88", label="TEM 130 ± 14 nm (lot)")
ax[0].axhline(TEM_D, color="0.4", lw=0.8)
xs = [1.0 / d for d in (800, 400, 200)]
for x, dl in zip(xs, (800, 400, 200)):
    ax[0].plot([x] * dls_lat[dl]["n"], [r["z_average"] for r in dls_lat[dl]["rs"]], "o", color=COL[0], alpha=0.7, ms=4)
ax[0].plot(xs, [dls_lat[d]["z"] for d in (800, 400, 200)], "-", color=COL[0], label="latex, 90° (3 repeats)")
ax[0].plot([0.01] * 3, [r["z_average"] for r in dls_lat[100]["rs"]], "^", color=COL[1], mfc="none", label="1:100, 173°")
ax[0].plot([0.01] * 3, [r["z_average"] for r in rs90], "s", color=COL[0], mfc="none", label="1:100, 90°")
ax[0].set_xscale("log")
ax[0].set_xlabel("Dilution (1/x)")
ax[0].set_ylabel("Z-average (nm)")
ax[0].set_ylim(110, 185)
ax[0].legend(fontsize=14, loc="lower center", bbox_to_anchor=(0.5, 1.02), ncol=2, frameon=False)
mark(ax[0], "a")
for dl, ls in zip((800, 400, 200), ("-", "--", ":")):
    a = m.dls_analyse(dls_lat[dl]["rs"][0]["rec"], m.DLS_MODELS[0])
    ax[1].semilogx(a["d"], a["p"] / a["p"].max(), ls, color=COL[0], label=f"1:{dl} (90°)")
a = m.dls_analyse(dls_lat[100]["rs"][0]["rec"], m.DLS_MODELS[0])
ax[1].semilogx(a["d"], a["p"] / a["p"].max(), "-", color=COL[1], label="1:100 (173°)")
ax[1].axvline(TEM_D, color="0.4", lw=0.8)
ax[1].set_xlim(40, 1500)
ax[1].set_xlabel("Diameter (nm)")
ax[1].set_ylabel("Intensity-weighted")
ax[1].legend(fontsize=14)
mark(ax[1], "b")
names = list(ratios)
ax[2].bar(range(3), [ratios[n] for n in names], color=COL[2], width=0.55)
ax[2].axhline(2.0, color="k", lw=0.9, ls="--")
ax[2].text(2.4, 2.03, "nominal 2", ha="right", fontsize=14)
ax[2].set_xticks(range(3))
ax[2].set_xticklabels(["UV A(532)", "DLS\nextinction", "DLS count\nrate / ND"], fontsize=16)
ax[2].set_ylabel("Ratio L2 / L1")
ax[2].set_ylim(0, 2.3)
mark(ax[2], "c")
axd["d"].bar(range(3), [R[800], R[400], R[200]], color=COL[3], width=0.55)
axd["d"].set_xticks(range(3))
axd["d"].set_xticklabels(["1:800", "1:400", "1:200"])
axd["d"].set_ylabel("Count rate / extinction")
axd["d"].set_xlabel("Dilution")
mark(axd["d"], "d")
fig.tight_layout(pad=1.6, h_pad=2.4)
fig.savefig(os.path.join(OUT, "fig2_latex_dls.png"), dpi=150)
plt.close(fig)

# =====================================================================================================================
# 3. Liposomes: UV against the simulator, with the DLS
# =====================================================================================================================
uv_lip = {}
for n, v in read_dsp("20251118 Lipo varie size").items():
    uv_lip["0.45" if "cuv" not in n else n.split(" in PBS ")[-1]] = v
LIP = {"100 nm": dict(key="0.1", dil=1.0, path=0.1, n=1), "200 nm": dict(key="0.2", dil=1.0, path=0.1, n=1),
       "450 nm (3 passes)": dict(key="0.45", dil=40.0, path=1.0, n=7)}
D18 = "20261118 lipo vuoti varie size DLS"
dls18 = read_dls(D18)
LIP_DLS = {"100 nm": "lipo 100 nm 1a10 in PBS", "200 nm": "lipo 200 nm 1a40 in PBS",
           "450 nm (3 passes)": ["lipo 450 nm 1a40 in PBS", "lipo 450 nm 1a40 in PBS bis"]}
lam_l = np.arange(350.0, 801.0, 8.0)


def lip_fit(key, dilution, path, n):
    ym = np.interp(lam_l, *uv_lip[key])

    def model(d):
        return m.vesicle_amplitude(lam_l, d, 0.2, n, 10.0 / dilution, path, acc_deg=10)

    def obj(d):
        a = model(d)
        k = float(a @ ym / (a @ a))
        return float(np.sqrt(np.mean((ym - k * a) ** 2))), k

    ds = np.geomspace(50, 1000, 12)
    r = [obj(d)[0] for d in ds]
    i = int(np.argmin(r))
    sol = minimize_scalar(lambda d: obj(d)[0], bounds=(ds[max(i - 1, 0)], ds[min(i + 1, len(ds) - 1)]),
                          method="bounded", options={"xatol": 1})
    rms, k = obj(sol.x)
    return float(sol.x), k, rms, model(sol.x), ym


lip = {}
fig, axd = plt.subplot_mosaic([["a"], ["b"], ["c"], ["d"]], figsize=(9, 21))
ax = [axd["a"], axd["b"], axd["c"]]
for i, (lab, p) in enumerate(LIP.items()):
    d, k, rms, A, ym = lip_fit(p["key"], p["dil"], p["path"], p["n"])
    dz = LIP_DLS[lab]
    zs = [r["z_average"] for key in (dz if isinstance(dz, list) else [dz]) for r in dls18[key]]
    pis = [r["pi"] for key in (dz if isinstance(dz, list) else [dz]) for r in dls18[key]]
    x, y = uv_lip[p["key"]]
    sel = (x >= 400) & (x <= 800)
    P = m.fit_power_law(x[sel], y[sel], m=0.0)[1]
    lip[lab] = dict(d_uv=d, k=k, rms=rms, lamellae=p["n"], z_dls=float(np.mean(zs)), pi_dls=float(np.mean(pis)), P=P)
    ax[i].plot(lam_l, ym, "-", color=COL[i], lw=1.6, label="measured")
    ax[i].plot(lam_l, k * A, "--", color="k", lw=1.0, label=f"simulator: {d:.0f} nm, {p['n']} lam., k = {k:.2f}")
    ax[i].set_title(f"{lab} (DLS {np.mean(zs):.0f} nm)")
    ax[i].set_xlabel("Wavelength (nm)")
    ax[i].legend(fontsize=14)
    mark(ax[i], "abc"[i])
ax[0].set_ylabel("Absorbance")
ax[2].set_ylabel("Absorbance")
x = np.arange(3)
axd["d"].bar(x - 0.2, [v["d_uv"] for v in lip.values()], 0.4, color=COL[2], label="UV simulator")
axd["d"].bar(x + 0.2, [v["z_dls"] for v in lip.values()], 0.4, color=COL[3], label="DLS Z-average")
axd["d"].set_xticks(x)
axd["d"].set_xticklabels(["100 nm", "200 nm", "450 nm\n(7 lamellae)"])
axd["d"].set_ylabel("Diameter (nm)")
axd["d"].legend(loc="upper left")
mark(axd["d"], "d")
fig.tight_layout(pad=1.6, h_pad=2.4)
fig.savefig(os.path.join(OUT, "fig3_liposome_uv_vs_simulator.png"), dpi=150)
plt.close(fig)
res["liposomes"] = lip

# =====================================================================================================================
# 4. Liposome '100' DLS against dilution, next to the latex
# =====================================================================================================================
lip100 = {}
for dl in (10, 20, 40):
    rs = dls19["lipo 100 1a%d in PBS" % dl]
    lip100[dl] = dict(z=float(np.mean([r["z_average"] for r in rs])), pi=float(np.mean([r["pi"] for r in rs])),
                      cr=float(np.mean([r["cr"] for r in rs])), rs=rs)
res["liposome_100_dls_19Nov"] = {dl: {x: v for x, v in d.items() if x != "rs"} for dl, d in lip100.items()}
fig, ax = plt.subplots(2, 1, figsize=(9, 11.5))
for dl in (10, 20, 40):
    ax[0].plot([1.0 / dl] * 3, [r["z_average"] for r in lip100[dl]["rs"]], "o", color=COL[3], alpha=0.7, ms=4)
ax[0].plot([1.0 / d for d in (40, 20, 10)], [lip100[d]["z"] for d in (40, 20, 10)], "-", color=COL[3], label="liposomes '100' (90°)")
ax[0].axhline(lip["100 nm"]["d_uv"], color="k", ls=":", lw=1.0, label="UV simulator, same lot (undiluted)")
ax[0].set_xscale("log")
ax[0].set_xlabel("Dilution (1/x)")
ax[0].set_ylabel("Z-average (nm)")
ax[0].legend(fontsize=14)
mark(ax[0], "a")
for dl, ls in zip((40, 20, 10), ("-", "--", ":")):
    a = m.dls_analyse(lip100[dl]["rs"][0]["rec"], m.DLS_MODELS[0])
    ax[1].semilogx(a["d"], a["p"] / a["p"].max(), ls, color=COL[3], label=f"1:{dl}")
ax[1].set_xlim(30, 2000)
ax[1].set_xlabel("Diameter (nm)")
ax[1].set_ylabel("Intensity-weighted")
ax[1].legend(fontsize=14)
mark(ax[1], "b")
fig.tight_layout(pad=1.6, h_pad=2.4)
fig.savefig(os.path.join(OUT, "fig4_liposome_dls.png"), dpi=150)
plt.close(fig)

# =====================================================================================================================
# 5. Theory: exponent P against diameter, with the measured exponents
# =====================================================================================================================
dgrid = np.geomspace(30, 800, 40)
lamP = np.arange(400.0, 801.0, 10.0)
P_solid = [m.fit_power_law(lamP, m.solid_sphere_amplitude(lamP, d, 0.1, 1.0, 1.0, n589=N589), m=0.0)[1] for d in dgrid]
P_ves = {s: m.vesicle_P_curve(dgrid, lamP, sigma=s) for s in (0.2, 0.3)}
res["P_theory"] = {"solid sphere (sigma 0.1)": dict(zip([round(float(d)) for d in dgrid], [round(float(p), 2) for p in P_solid]))}
fig, ax = plt.subplots(figsize=(9, 7.5))
ax.semilogx(dgrid, P_solid, "-", color=COL[0], label="solid polystyrene spheres (σ = 0.1)")
for s, ls in zip((0.2, 0.3), ("-", "--")):
    ax.semilogx(dgrid, P_ves[s], ls, color=COL[3], label=f"unilamellar vesicles (σ = {s})")
for k in dil:
    ax.plot(TEM_D, lat[k]["P"], "o", color=COL[0], ms=5)
ax.plot([], [], "o", color=COL[0], label="latex L1-L3 measured P at the TEM diameter")
for lab, v in lip.items():
    ax.plot(v["z_dls"], v["P"], "s", color=COL[3], ms=6)
ax.plot([], [], "s", color=COL[3], label="liposomes measured P at the DLS Z-average")
ax.set_xlabel("Diameter (nm)")
ax.set_ylabel("Exponent P (400-800 nm)")
ax.set_xlim(30, 800)
ax.set_xticks([30, 50, 100, 200, 400, 800])
ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
ax.minorticks_off()
ax.legend(fontsize=14, loc="lower left")
fig.tight_layout(pad=1.6, h_pad=2.4)
fig.savefig(os.path.join(OUT, "fig5_P_vs_size.png"), dpi=150)
plt.close(fig)

# ---- print the numbers --------------------------------------------------------------------------------------------
import json
print(json.dumps(res, indent=1, default=lambda o: float(o) if isinstance(o, (np.floating, np.integer)) else str(o)))
