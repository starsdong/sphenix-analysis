#!/usr/bin/env python3
"""make_numbers.py <kappa csv> [287|286] -- all numbers, tables and PASS/FAIL checks of the kappa-calibration note.
Reads logs/ (copied/produced by build.sh) and writes tables/numbers.tex (\\newcommand, letters-only names) and tables/tab_*.tex."""
import csv, os, re, sys

NOTE = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
assert os.path.exists(os.path.join(NOTE, "main.tex")) and os.path.exists(os.path.join(NOTE, "build.sh"))
CSV = sys.argv[1]
MODE = int(sys.argv[2]) if len(sys.argv) > 2 else 287
LAST = (81601, 5030)   # last segment; excluded in MODE 286 (fits not yet redone with it)
LOG = lambda f: open(os.path.join(NOTE, "logs", f)).read()
mac = {}


def m(n, v):
    assert re.fullmatch(r"[A-Za-z]+", n), n
    assert n not in mac, f"duplicate macro {n}"
    mac[n] = str(v)


def num(x):
    return f"{x:,}".replace(",", "{,}")


PDG = {"Ks": 497.611, "Lambda": 1115.683}
PER = {"A": "PerA", "B": "PerB", "all": "PerAll"}
VER = {"kfp": "Kfp", "raw": "Raw", "kappa": "Kap"}

# ---------------------------------------------------------------- kappa table
rows = [l for l in open(CSV) if not l.startswith("#")]
r = list(csv.DictReader(rows))
nodes = {k: len(set(x[k] for x in r)) for k in ("q", "pT", "eta", "phi")}
m("TabRows", num(len(r))); m("TabNq", nodes["q"]); m("TabNpt", nodes["pT"]); m("TabNeta", nodes["eta"]); m("TabNphi", nodes["phi"])
pts = sorted(set(float(x["pT"]) for x in r))
m("TabPtMin", f"{pts[0]:g}"); m("TabPtMax", f"{pts[-1]:g}")
etas = sorted(set(float(x["eta"]) for x in r)); m("TabEtaMin", f"{etas[0]:g}"); m("TabEtaMax", f"{etas[-1]:g}")
d = {(float(x["q"]), float(x["pT"]), float(x["eta"]), round(float(x["phi"]), 4)): float(x["kappa"]) for x in r}
for q, w in ((-1.0, "Minus"), (1.0, "Plus")):
    for p, pw in ((0.3, "PtA"), (1.0, "PtB"), (2.0, "PtC")):
        m(f"TabK{w}{pw}", f"{d[(q, p, 0.0, 0.0)]:.4f}")
kap = [float(x["kappa"]) for x in r]
m("TabKMin", f"{min(kap):.3f}"); m("TabKMax", f"{max(kap):.3f}")
head = "\n".join(l for l in open(CSV) if l.startswith("#"))
m("TabXref", re.search(r"x_ref = (\S+) ln", head).group(1)); m("TabXscale", re.search(r"x_scale = ([\d.]+)", head).group(1))
m("TabXclamp", re.search(r"x_clamp = ([\d.]+)", head).group(1))
# kappa vs phi spread at pT 0.5, eta 0 (illustration of the phi dependence)
for q, w in ((-1.0, "Minus"), (1.0, "Plus")):
    v = [d[(q, 0.5, 0.0, ph)] for (qq, pp, ee, ph) in d if qq == q and pp == 0.5 and ee == 0.0]
    m(f"TabPhiLo{w}", f"{min(v):.3f}"); m(f"TabPhiHi{w}", f"{max(v):.3f}")
    v = [d[(q, 0.5, ee, 0.0)] for (qq, pp, ee, ph) in d if qq == q and pp == 0.5 and ph == 0.0]
    m(f"TabEtaLo{w}", f"{min(v):.3f}"); m(f"TabEtaHi{w}", f"{max(v):.3f}")
# value table at eta = 0, phi = 0
vt = []
for p in (0.1, 0.2, 0.3, 0.5, 0.7, 1.0, 1.5, 2.0, 3.0, 5.0):
    vt.append(f"{p:g} & {d[(-1.0, p, 0.0, 0.0)]:.4f} & {d[(1.0, p, 0.0, 0.0)]:.4f} \\\\")

# ---------------------------------------------------------------- validation logs
cmpf = LOG("cmpKappaFill.txt")
nh = len(re.findall(r"^MAXDIFF", cmpf, re.M)); nzero = len(re.findall(r"^MAXDIFF \S+\s+\S+\s+0 ", cmpf, re.M))
m("ChkFillN", nh); m("ChkFill", "PASS" if "CHECK PASS" in cmpf and nh == nzero == 18 else "FAIL")
rp = re.search(r"REPRO compared (\d+) identical (\d+) differ (\d+) no-old (\d+)", LOG("repro_287.txt"))
m("ReproN", rp.group(1)); m("ReproSame", rp.group(2))
m("ChkRepro", "PASS" if rp.group(1) == rp.group(2) == "287" else "FAIL")
fl = LOG("fill_validate.txt")
sp = re.search(r"kappa\(q=-1,pT=1,eta=0,phi=0\) = (\S+), kappa\(q=\+1,pT=1,eta=0,phi=0\) = (\S+)", fl)
miss = int(re.search(r"missing (\d+)", fl).group(1))
m("ChkTable", "PASS" if miss == 0 and len(r) == nodes["q"] * nodes["pT"] * nodes["eta"] * nodes["phi"] else "FAIL")
m("ChkInterp", "PASS" if abs(float(sp.group(1)) - d[(-1.0, 1.0, 0.0, 0.0)]) < 1e-5 and abs(float(sp.group(2)) - d[(1.0, 1.0, 0.0, 0.0)]) < 1e-5 else "FAIL")

# ---------------------------------------------------------------- jobs, runs, periods
jobs = []
for l in LOG("joblines.txt").splitlines():
    g = re.match(r"SEG (\d+)-(\d+) Ks QA candidates (\d+), Lambda QA candidates (\d+), events (\d+) \| mean kappa \(QA\): "
                 r"Ks d1 (\S+) d2 (\S+) ; Lambda d1 (\S+) d2 (\S+) \| (\d+)", l)
    assert g, l
    jobs.append(dict(run=int(g.group(1)), seg=int(g.group(2)), nks=int(g.group(3)), nla=int(g.group(4)), nev=int(g.group(5)),
                     k=[float(g.group(i)) for i in range(6, 10)], ok=int(g.group(10))))
m("JobsN", len(jobs)); m("JobsOk", sum(j["ok"] for j in jobs))
m("ChkJobs", "PASS" if len(jobs) == 287 and all(j["ok"] == 1 for j in jobs) else "FAIL")
if MODE == 286:
    jobs = [j for j in jobs if (j["run"], j["seg"]) != LAST]
m("NSegUsed", len(jobs))
inA = lambda run: 79509 <= run <= 79516
runs = sorted(set(j["run"] for j in jobs))
m("NRuns", len(runs)); m("NRunsA", sum(inA(x) for x in runs)); m("NRunsB", sum(not inA(x) for x in runs))
m("RunFirst", runs[0]); m("RunLast", runs[-1])
runsB = [x for x in runs if not inA(x)]
m("RunBFirst", runsB[0]); m("RunBLast", runsB[-1])
rrows = []
for p, sel in (("A", inA), ("B", lambda x: not inA(x)), ("all", lambda x: True)):
    J = [j for j in jobs if sel(j["run"])]
    P = PER[p]
    m(f"{P}Seg", len(J)); m(f"{P}Ev", num(sum(j["nev"] for j in J)))
    m(f"{P}NKs", num(sum(j["nks"] for j in J))); m(f"{P}NLa", num(sum(j["nla"] for j in J)))
    nks = sum(j["nks"] for j in J)
    for i, w in enumerate(("KsOne", "KsTwo", "LaOne", "LaTwo")):
        wt = [j["nks"] if i < 2 else j["nla"] for j in J]
        m(f"{P}Kap{w}", f"{sum(j['k'][i] * x for j, x in zip(J, wt)) / sum(wt):.4f}")
for run in runs:
    J = [j for j in jobs if j["run"] == run]
    nev = sum(j["nev"] for j in J); nks = sum(j["nks"] for j in J); nla = sum(j["nla"] for j in J)
    k = [sum(j["k"][i] * (j["nks"] if i < 2 else j["nla"]) for j in J) / max(nks if i < 2 else nla, 1) for i in range(4)]
    segs = sorted(j["seg"] for j in J)
    rrows.append(f"{run} & {'A' if inA(run) else 'B'} & {len(J)} & {num(nev)} & {num(nks)} & {nks/nev:.3f} & {num(nla)} & "
                 f"{k[0]:.4f} & {k[1]:.4f} & {k[2]:.4f} & {k[3]:.4f} \\\\")
segA = sorted(j["seg"] for j in jobs if j["run"] == 79509)
m("SegsScanAnine", ", ".join(str(s) for s in segA if s > 30))
ksev = {}
for run in runs:
    J = [j for j in jobs if j["run"] == run]
    ksev[run] = sum(j["nks"] for j in J) / sum(j["nev"] for j in J)
for p, sel in (("A", inA), ("B", lambda x: not inA(x))):
    v = [ksev[x] for x in runs if sel(x)]
    m(f"{PER[p]}KsEvMin", f"{min(v):.3f}"); m(f"{PER[p]}KsEvMax", f"{max(v):.3f}")

# ---------------------------------------------------------------- fits
def fits(tag):
    t = LOG(f"fitMomRes_{tag}.txt")
    out = {}
    for l in t.splitlines():
        if l.startswith("FIT "):
            w = l.split()
            out[w[1]] = dict(a=float(w[3]), alo=float(w[4]), ahi=float(w[5]), b=float(w[7]), blo=float(w[8]), bhi=float(w[9]),
                             s=float(w[11]), slo=float(w[12]), shi=float(w[13]), chi2=float(w[15]), ndf=int(w[17]),
                             edge=int(w[19]) + int(w[21]) + int(w[23]))
    cr = re.search(r"CROSS chi2\(Lambda \| Ks best\) (\S+) / (\d+) ; chi2\(Lambda \| own best\) (\S+)", t)
    out["dchiLa"] = float(cr.group(1)) - float(cr.group(3))
    out["bins"] = []
    for l in t.splitlines():
        if l.startswith("BIN "):
            w = l.split()
            out["bins"].append(dict(ch=w[1], lo=float(w[2]), hi=float(w[3]), chi2=float(w[5]), ndf=int(w[7]), N=float(w[9]),
                                    gd=float(w[11]), gde=float(w[12]), gs=float(w[14])))
    return out


F = {}
PERIODS = ("A", "B", "all") if MODE == 287 else ("A", "B")
for p in PERIODS:
    for v in ("kfp", "raw", "kappa"):
        F[(p, v)] = fits(f"k287_{p}_{v}")
        f = F[(p, v)]
        P = PER[p] + VER[v]
        for ch, cw in (("Ks", "Ks"), ("Lambda", "La"), ("joint", "Jt")):
            x = f[ch]
            m(f"{P}{cw}A", f"{100*x['a']:.2f}"); m(f"{P}{cw}AUp", f"{100*(x['ahi']-x['a']):.2f}"); m(f"{P}{cw}ADn", f"{100*(x['a']-x['alo']):.2f}")
            m(f"{P}{cw}B", f"{100*x['b']:.2f}"); m(f"{P}{cw}S", f"{x['s']:.3f}")
            m(f"{P}{cw}Chi", f"{x['chi2']:.1f}"); m(f"{P}{cw}Ndf", x["ndf"])
        m(f"{P}DchiLa", f"{f['dchiLa']:.1f}")
edges = [(k, c) for k, f in F.items() for c in ("Ks", "joint") if f[c]["edge"]]
m("ChkEdge", "PASS" if not edges else "FAIL")
laEdges = [f"{p}/{v}" for (p, v), f in F.items() if f["Lambda"]["edge"]]
m("LaEdgeList", ", ".join(laEdges) if laEdges else "none")
lab = {"kfp": "KFP mass", "raw": "recomputed, no $\\kappa$", "kappa": "recomputed, with $\\kappa$"}
pm = lambda v, lo, hi, sc, nd: f"${sc*v:.{nd}f}^{{+{sc*(hi-v):.{nd}f}}}_{{-{sc*(v-lo):.{nd}f}}}$"
frows = []
for p in PERIODS:
    for v in ("kfp", "raw", "kappa"):
        f = F[(p, v)]
        for ch, cl in (("Ks", "\\Ks"), ("Lambda", "\\La"), ("joint", "joint")):
            x = f[ch]
            first = f"{p} & {lab[v]}" if ch == "Ks" else " & "
            dch = f"{f['dchiLa']:.1f}" if ch == "Ks" else ""
            frows.append(f"{first} & {cl} & {pm(x['a'], x['alo'], x['ahi'], 100, 2)} & {pm(x['b'], x['blo'], x['bhi'], 100, 2)} & "
                         f"{pm(x['s'], x['slo'], x['shi'], 1, 4)} & {x['chi2']:.1f}/{x['ndf']} & {dch} & {'yes' if x['edge'] else 'no'} \\\\")
        frows.append("\\addlinespace")
    frows[-1] = "\\midrule"
frows.pop()
# compact summary table: Ks (= joint) per period and version
srows = []
for p in PERIODS:
    for v in ("kfp", "raw", "kappa"):
        x = F[(p, v)]["Ks"]; j = F[(p, v)]["joint"]
        same = abs(x["a"] - j["a"]) < 1e-9 and abs(x["b"] - j["b"]) < 1e-9 and abs(x["s"] - j["s"]) < 1e-9
        srows.append(f"{p if v == 'kfp' else ''} & {lab[v]} & {100*x['a']:.2f} & {100*x['b']:.2f} & {x['s']:.3f} & {x['chi2']:.1f}/{x['ndf']} & "
                     f"{F[(p, v)]['dchiLa']:.1f} & {'yes' if same else 'no'} \\\\")
    srows.append("\\midrule")
srows.pop()
# per-bin data statistics and Gaussian widths, kappa version (template fit at the Ks best point)
brows = []
for p in ("A", "B"):
    for b in F[(p, "kappa")]["bins"]:
        braw = [x for x in F[(p, "raw")]["bins"] if x["ch"] == b["ch"] and x["lo"] == b["lo"]][0]
        cl = "\\Ks" if b["ch"] == "Ks" else "\\La"
        brows.append(f"{p} & {cl} & {b['lo']:.1f}--{b['hi']:.1f} & {braw['chi2']:.1f}/{braw['ndf']} & {b['chi2']:.1f}/{b['ndf']} \\\\")
    brows.append("\\midrule")
brows.pop()
# change in a and s
for p in PERIODS:
    m(f"{PER[p]}ARatio", f"{F[(p,'kappa')]['Ks']['a']/F[(p,'raw')]['Ks']['a']:.2f}")
# template statistics
fl = open(os.path.join(NOTE, "logs", "fit_k287_A_kappa.log")).read()
m("TmplKs", num(int(re.search(r"Ks: (\d+) DDSim candidates", fl).group(1))))
m("TmplLa", num(int(re.search(r"Lambda: (\d+) DDSim candidates", fl).group(1))))
ent = {}
for p in PERIODS:
    for v in ("kfp", "raw", "kappa"):
        t = open(os.path.join(NOTE, "logs", f"fit_k287_{p}_{v}.log")).read()
        ent[(p, v, "Ks")] = [int(x) for x in re.search(r"Ks data entries per pT bin \(fit range\): ([\d ]+)", t).group(1).split()]
        ent[(p, v, "La")] = [int(x) for x in re.search(r"Lambda data entries per pT bin \(fit range\): ([\d ]+)", t).group(1).split()]
m("AllFitKs", num(sum(ent[('A', 'kfp', 'Ks')]) + sum(ent[('B', 'kfp', 'Ks')]))); m("AllFitLa", num(sum(ent[('A', 'kfp', 'La')]) + sum(ent[('B', 'kfp', 'La')])))
erows = []
ptK = ["0.2--0.5", "0.5--0.7", "0.7--0.9", "0.9--1.2", "1.2--1.6", "1.6--2.2", "2.2--4.0"]
ptL = ["0.5--0.9", "0.9--1.2", "1.2--1.6", "1.6--2.4", "2.4--4.0"]
for ch, pts_, cl in (("Ks", ptK, "\\Ks"), ("La", ptL, "\\La+$\\bar\\Lambda$")):
    for i, b in enumerate(pts_):
        erows.append(f"{cl if i == 0 else ''} & {b} & " + " & ".join(num(ent[(p, v, ch)][i]) for p in ("A", "B") for v in ("kfp", "kappa")) + " \\\\")
    erows.append("\\midrule")
erows.pop()

# ---------------------------------------------------------------- Gaussian peaks (peaks.C)
pk = {}
for g in re.finditer(r"^GPEAK (\S+) (\S+) (\S+) (\S+) (\S+) (\S+) (\S+) (\S+) (\S+) (\S+) (\d)", LOG("peaks.txt"), re.M):
    pk[(g.group(1), g.group(2), g.group(3), float(g.group(4)), float(g.group(5)))] = (
        float(g.group(6)), float(g.group(7)), float(g.group(8)), float(g.group(9)), int(g.group(11)))
prow = {"Ks": [], "Lambda": []}
for ch, full in (("Ks", (0.2, 4.0)), ("Lambda", (0.5, 4.0))):
    bins = sorted({(k[3], k[4]) for k in pk if k[0] == ch})
    bins = [b for b in bins if b != full] + [full]
    for b in bins:
        cells = []
        for p in ("A", "B"):
            for v in ("raw", "kappa"):
                x = pk[(ch, p, v) + b]
                cells += [f"{x[0]:.2f}$\\pm${x[1]:.2f}", f"{x[2]:.2f}$\\pm${x[3]:.2f}"] if x[4] else ["--", "--"]
        prow[ch].append(f"{b[0]:.1f}--{b[1]:.1f}{' (all)' if b == full else ''} & " + " & ".join(cells) + " \\\\")
    c = "Ks" if ch == "Ks" else "La"
    for p in ("A", "B"):
        for v in ("raw", "kappa"):
            x = pk[(ch, p, v) + full]
            P = PER[p] + VER[v]
            m(f"{P}{c}GM", f"{x[0]:.2f}"); m(f"{P}{c}GMe", f"{x[1]:.2f}"); m(f"{P}{c}GS", f"{x[2]:.2f}"); m(f"{P}{c}GSe", f"{x[3]:.2f}")
            m(f"{P}{c}GPull", f"{(x[0]-PDG[ch])/x[1]:+.1f}")
            m(f"{P}{c}Below", f"{PDG[ch]-x[0]:.2f}")
    m(f"Pdg{c}", f"{PDG[ch]:.3f}")
# Ks width at 1.6-2.2
for p in ("A", "B"):
    for v in ("raw", "kappa"):
        x = pk[("Ks", p, v, 1.6, 2.2)]
        m(f"{PER[p]}{VER[v]}KsWidHigh", f"{x[2]:.1f}")
# ---------------------------------------------------------------- mass distributions (massdist.C), pT-integrated
mdn = {}
for g in re.finditer(r"^MDIST (\S+) (\S+) (\S+) (\S+) (\S+) (\S+) (\S+)", LOG("massdist.txt"), re.M):
    if (g.group(1) == "Ks" and g.group(4) == "0.2" and g.group(5) == "4.0") or (g.group(1) == "Lambda" and g.group(4) == "0.5" and g.group(5) == "4.0"):
        mdn[(g.group(1), g.group(2), g.group(3))] = float(g.group(6))
for ch, c in (("Ks", "Ks"), ("Lambda", "La")):
    for p in ("A", "B"):
        P = PER[p]
        m(f"{P}{c}MdKfp", num(int(mdn[(ch, p, 'kfp')]))); m(f"{P}{c}MdRaw", num(int(mdn[(ch, p, 'raw')])))
        m(f"{P}{c}MdKap", num(int(mdn[(ch, p, 'kappa')])))
        m(f"{P}{c}MdRatio", f"{mdn[(ch, p, 'raw')] / mdn[(ch, p, 'kfp')]:.2f}")
# ---------------------------------------------------------------- resolution curves (plot_res.C)
rs = LOG("res.txt")
for g in re.finditer(r"^CURVE k287_(\w+)_(\w+) (\w+) (\S+) (\S+) (\S+)", rs, re.M):
    P = PER[g.group(1)] + VER[g.group(2)]
    s_ = {"pi": "Pi", "p": "P"}[g.group(3)]
    for i, w in ((4, "One"), (5, "Two"), (6, "Three")):
        m(f"{P}Res{s_}{w}", f"{100*float(g.group(i)):.1f}")
single = []
for g in re.finditer(r"^BAND k287_(\w+?)_(\w+) pi 2.0 (\S+) (\S+)", rs, re.M):
    if float(g.group(3)) == float(g.group(4)):
        single.append(f"{g.group(1)} {'with' if g.group(2) == 'kappa' else 'no'} $\\kappa$")
m("ResBandNote", ("only the best node passes for " + ", ".join(single) + ", so those bands are hidden under their lines")
  if single else "every fit has more than one node inside the band")
# ---------------------------------------------------------------- run-group fits (KFP mass)
gl = []
for tag, label in (("kfp31", "79509 seg.\\ 0--30"), ("k287_A_kfp", "79509--79516 (A)"), ("g7951x", "79510--79516"),
                   ("g809xx", "80904--80915"), ("g81xxx", "81556--81609"), ("k287_B_kfp", "B (all non-A)")) + ((("k287_all_kfp", "all"),) if MODE == 287 else ()):
    f = fits(tag)
    x = f["Ks"]
    gl.append(f"{label} & {pm(x['a'], x['alo'], x['ahi'], 100, 2)} & {pm(x['b'], x['blo'], x['bhi'], 100, 2)} & {x['s']:.3f} & {x['chi2']:.1f}/{x['ndf']} \\\\")
m("GrpNine", f"{100*fits('g809xx')['Ks']['a']:.2f}"); m("GrpEightOne", f"{100*fits('g81xxx')['Ks']['a']:.2f}")
m("GrpThirtyOne", f"{100*fits('kfp31')['Ks']['a']:.2f}")
# ---------------------------------------------------------------- 31-seg kappa profile (local, from the earlier study)
kp = {}
for g in re.finditer(r"^KAPPA (\S+) (\S+) (\S+) (\S+)", LOG("diagKappa.txt"), re.M):
    kp[(g.group(1), float(g.group(2)))] = (float(g.group(3)), int(float(g.group(4))))
ptv = sorted(set(k[1] for k in kp))
krows = []
for pn, l in (("pKappa_Ks_d1", "\\Ks: $\\pi^-$"), ("pKappa_Ks_d2", "\\Ks: $\\pi^+$"), ("pKappa_La_d1", "\\La+$\\bar\\Lambda$: $\\pi$"),
              ("pKappa_La_d2", "\\La+$\\bar\\Lambda$: p, $\\bar{\\rm p}$")):
    krows.append(f"{l} & " + " & ".join(f"{kp[(pn, p)][0]:.3f}" if kp[(pn, p)][1] >= 50 else "--" for p in ptv) + " \\\\")
for ch, h in (("Ks", "hDm_Ks"), ("La", "hDm_La")):
    g = re.search(rf"^DM {h} (\S+) (\S+)", LOG("diagKappa.txt"), re.M)
    m(f"Dm{ch}Mean", f"{float(g.group(1)):+.2f}"); m(f"Dm{ch}Rms", f"{float(g.group(2)):.1f}")
m("KprofPts", " & ".join(f"{p:g}" for p in ptv))
if MODE == 286:
    m("SampleStatus", "Status: period B and the combined fits use 191 of the 192 period-B segments (286 of 287); the last segment "
      "(run 81601, segment 5030) finished after the fits were started, and the combined A+B fits are still running. This version "
      "will be updated.")
else:
    m("SampleStatus", "All 287 segments are used.")

# ---------------------------------------------------------------- tight Lambda selection (Sec. tight)
sc = {}
for l in LOG("scanLambdaCuts.txt").splitlines():
    if l.startswith("SCAN "):
        a_, b_, c_ = l.split("|")
        k = tuple(a_.split()[1:])
        pbv = list(map(float, b_.split())); S_, B_, sb_, se_ = map(float, c_.split())
        sc[k] = dict(S=S_, B=B_, sb=sb_, se=se_, pb=pbv)
m("ScanNCuts", len(sc))
for key, w in ((("0.02", "0.02", "0.5", "0.99500"), "Base"), (("0.10", "0.03", "0.5", "0.99950"), "Sel"),
               (("0.15", "0.04", "0.5", "0.99970"), "Alt"), (("0.02", "0.02", "0.5", "0.99950"), "Dira"),
               (("0.02", "0.02", "2.0", "0.99500"), "Lxy")):
    x = sc[key]
    m(f"Scan{w}S", f"{x['S']:.0f}"); m(f"Scan{w}B", f"{x['B']:.0f}"); m(f"Scan{w}SB", f"{x['sb']:.2f}"); m(f"Scan{w}Seff", f"{x['se']:.0f}")
    m(f"Scan{w}Keep", f"{100*x['S']/sc[('0.02','0.02','0.5','0.99500')]['S']:.0f}")
m("ScanGain", f"{sc[('0.10','0.03','0.5','0.99950')]['se']/sc[('0.02','0.02','0.5','0.99500')]['se']:.2f}")
best = max(sc.items(), key=lambda kv: kv[1]["se"])
m("ScanBestCut", f"DIRA$>${float(best[0][3]):g}, $|{{\\rm DCA}}_{{xy}}(\\pi)|>{float(best[0][0]):g}$, $|{{\\rm DCA}}_{{xy}}(p)|>{float(best[0][1]):g}$~cm")
m("ScanBestSeff", f"{best[1]['se']:.0f}")
bt = LOG("cutBiasTest.txt")
for ch, c in (("Ks", "Ks"), ("Lambda", "La")):
    g = re.search(rf"^BIAS {ch} pt (\S+)-4.0 data: sigma base (\S+)\+-(\S+) tight (\S+)\+-(\S+) ratio (\S+)\+-(\S+) \(kept (\S+) of entries\) \| DDSim: base (\S+) tight (\S+) ratio (\S+)", bt, re.M)
    if ch == "Ks":   # integrated row is the one starting at 0.2 that ends at 4.0 and is last; take the last match
        g = list(re.finditer(rf"^BIAS {ch} pt (\S+)-4.0 data: sigma base (\S+)\+-(\S+) tight (\S+)\+-(\S+) ratio (\S+)\+-(\S+) \(kept (\S+) of entries\) \| DDSim: base (\S+) tight (\S+) ratio (\S+)", bt, re.M))[-1]
    else:
        g = list(re.finditer(rf"^BIAS {ch} pt (\S+)-4.0 data: sigma base (\S+)\+-(\S+) tight (\S+)\+-(\S+) ratio (\S+)\+-(\S+) \(kept (\S+) of entries\) \| DDSim: base (\S+) tight (\S+) ratio (\S+)", bt, re.M))[-1]
    m(f"Bias{c}RatD", g.group(6)); m(f"Bias{c}RatDe", g.group(7)); m(f"Bias{c}RatS", g.group(11))
biasrows = []
for g in re.finditer(r"^BIAS (\S+) pt (\S+)-(\S+) data: sigma base (\S+)\+-(\S+) tight (\S+)\+-(\S+) ratio (\S+)\+-(\S+) \(kept (\S+) of entries\) \| DDSim: base (\S+) tight (\S+) ratio (\S+)", bt, re.M):
    biasrows.append(f"{'\\Ks' if g.group(1) == 'Ks' else '\\La+$\\bar\\Lambda$'} & {g.group(2)}--{g.group(3)} & {g.group(4)}$\\pm${g.group(5)} & {g.group(6)}$\\pm${g.group(7)} & "
                    f"{g.group(8)}$\\pm${g.group(9)} & {g.group(11)} & {g.group(12)} & {g.group(13)} \\\\")
rv = re.search(r"REPROV2 segments (\d+) differing-histograms (\d+) missing-old (\d+)", LOG("repro_v2.txt"))
m("ReproTwoN", rv.group(1)); m("ChkReproTwo", "PASS" if rv.groups() == ("287", "0", "0") else "FAIL")
ntl = int(re.search(r"Lambda: (\d+) DDSim candidates \(passAll \+ tight", LOG("fit_k287t_A_kappa.log")).group(1))
m("TmplLaTight", num(ntl)); m("TmplLaTightFrac", f"{ntl / int(re.search(r'Lambda: (\d+) DDSim candidates', fl).group(1)):.2f}")
# Gaussian peaks QA vs tight
lp = {}
for g in re.finditer(r"^LPEAK (\S+) (\S+) (\S+) (\S+) (\S+) mean (\S+) (\S+) sigma (\S+) (\S+) S (\S+) B2s (\S+) SB (\S+)", LOG("laTightPeaks.txt"), re.M):
    lp[(g.group(1), g.group(2), g.group(3), float(g.group(4)), float(g.group(5)))] = [float(g.group(i)) for i in range(6, 13)]
lrows = []
for p in ("A", "B", "all"):
    for v in ("raw", "kappa"):
        q, t = lp[(p, v, "qa", 0.5, 4.0)], lp[(p, v, "tight", 0.5, 4.0)]
        lrows.append(f"{p if v == 'raw' else ''} & {'no $\\kappa$' if v == 'raw' else 'with $\\kappa$'} & {q[0]:.2f}$\\pm${q[1]:.2f} & {t[0]:.2f}$\\pm${t[1]:.2f} & "
                     f"{q[2]:.2f} & {t[2]:.2f} & {q[6]:.2f} & {t[6]:.2f} & {100*t[4]/q[4]:.0f} \\\\")
        P = PER[p] + VER[v]
        m(f"{P}LpQaM", f"{q[0]:.2f}"); m(f"{P}LpQaMe", f"{q[1]:.3f}"); m(f"{P}LpTiM", f"{t[0]:.2f}"); m(f"{P}LpTiMe", f"{t[1]:.3f}")
        m(f"{P}LpQaSB", f"{q[6]:.2f}"); m(f"{P}LpTiSB", f"{t[6]:.2f}"); m(f"{P}LpKeep", f"{100*t[4]/q[4]:.0f}")
        m(f"{P}LpTiPull", f"{(t[0]-PDG['Lambda'])/t[1]:+.1f}"); m(f"{P}LpTiAbove", f"{t[0]-PDG['Lambda']:.2f}")
    lrows.append("\\midrule")
lrows.pop()
# template fits QA vs tight (Lambda-only, joint, Delta chi2)
T = {}
trows = []
for p in ("A", "B", "all"):
    for v in ("raw", "kappa"):
        T[(p, v)] = fits(f"k287t_{p}_{v}")
        for cut, f in (("QA", F[(p, v)]), ("tight", T[(p, v)])):
            x = f["Lambda"]; j = f["joint"]
            trows.append(f"{p if (v == 'raw' and cut == 'QA') else ''} & {('no $\\kappa$' if v == 'raw' else 'with $\\kappa$') if cut == 'QA' else ''} & {cut} & "
                         f"{pm(x['a'], x['alo'], x['ahi'], 100, 2)} & {100*x['b']:.2f} & {pm(x['s'], x['slo'], x['shi'], 1, 4)} & {x['chi2']:.1f}/{x['ndf']} & "
                         f"{f['dchiLa']:.1f} & {100*j['a']:.2f}, {100*j['b']:.2f}, {j['s']:.3f} \\\\")
            P = PER[p] + VER[v] + ("Ti" if cut == "tight" else "Qa")
            m(f"{P}LaA", f"{100*x['a']:.2f}"); m(f"{P}LaS", f"{x['s']:.3f}"); m(f"{P}LaSLo", f"{x['slo']:.4f}"); m(f"{P}LaSHi", f"{x['shi']:.4f}")
            m(f"{P}LaSW", f"{(x['shi']-x['slo'])/2*1000:.1f}"); m(f"{P}LaChi", f"{x['chi2']:.1f}"); m(f"{P}Dchi", f"{f['dchiLa']:.1f}")
        trows.append("\\addlinespace")
    trows[-1] = "\\midrule"
trows.pop()
tedge = [k for k, f in T.items() for c in ("Ks", "Lambda", "joint") if f[c]["edge"]]
m("ChkEdgeTight", "PASS" if not tedge else "FAIL")
same = all(abs(T[k]["Ks"][q] - F[k]["Ks"][q]) < 1e-9 for k in T for q in ("a", "b", "s"))
m("ChkKsSame", "PASS" if same else "FAIL")
# per-bin Gaussian widths data vs template at the Ks best point, B with kappa, QA vs tight
wrows = []
for p in ("A", "B"):
    for cut, f in (("QA", F[(p, "kappa")]), ("tight", T[(p, "kappa")])):
        bl = [b for b in f["bins"] if b["ch"] == "Lambda"]
        wrows.append(f"{p if cut == 'QA' else ''} & {cut} & " + " & ".join(
            (f"{b['gd']:.2f}$\\pm${b['gde']:.2f} / {b['gs']:.2f}" if (b["gde"] < 1.0 and b["gd"] > 0.35) else f"-- / {b['gs']:.2f}") for b in bl) + " \\\\")
    wrows.append("\\midrule")
wrows.pop()
tb = [b for b in T[("B", "kappa")]["bins"] if b["ch"] == "Lambda"]; qb = [b for b in F[("B", "kappa")]["bins"] if b["ch"] == "Lambda"]
for i, w in ((2, "One"), (3, "Two")):
    m(f"WidTi{w}", f"{(tb[i]['gd'] - tb[i]['gs']) / tb[i]['gde']:.1f}"); m(f"WidQa{w}", f"{(qb[i]['gd'] - qb[i]['gs']) / qb[i]['gde']:+.1f}")
# ---------------------------------------------------------------- write
W = lambda f, s: open(os.path.join(NOTE, "tables", f), "w").write("% generated by scripts/make_numbers.py\n" + s + "\n")
W("tab_kappa_nodes.tex", "\n".join(vt))
W("tab_runs.tex", "\n".join(rrows))
W("tab_fits.tex", "\n".join(frows))
W("tab_summary.tex", "\n".join(srows))
W("tab_binchi.tex", "\n".join(brows))
W("tab_entries.tex", "\n".join(erows))
W("tab_peaks_ks.tex", "\n".join(prow["Ks"]))
W("tab_peaks_la.tex", "\n".join(prow["Lambda"]))
W("tab_groups.tex", "\n".join(gl))
W("tab_kprof.tex", "\n".join(krows))
W("tab_tight_bias.tex", "\n".join(biasrows))
W("tab_tight_peaks.tex", "\n".join(lrows))
W("tab_tight_fits.tex", "\n".join(trows))
W("tab_tight_widths.tex", "\n".join(wrows))
W("numbers.tex", "\n".join(f"\\newcommand{{\\{k}}}{{{v}}}" for k, v in mac.items()))
for k, v in mac.items():
    print(k, v)
print("CHECKS", {k: v for k, v in mac.items() if k.startswith("Chk")})
