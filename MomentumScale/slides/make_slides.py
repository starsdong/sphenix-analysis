#!/usr/bin/env python3
"""kappa / tight-Lambda momentum-calibration summary deck in Xin Dong's slide style (skill my-style-pptx).
All numbers come from the note's script output (../logs/numbers.txt, written by scripts/make_numbers.py).
run: env -u PYTHONPATH ~/.hermes/cache/scratch/pptxenv/bin/python make_slides.py
"""
import os, re, sys
sys.path.insert(0, os.path.expanduser('~/.hermes/skills/productivity/my-style-pptx/scripts'))
from mystyle import Deck, RED, BLUE, GREY
from pptx.util import Inches, Pt
from pptx.dml.color import RGBColor
from pptx.enum.text import PP_ALIGN

HERE = os.path.dirname(os.path.abspath(__file__))
NOTE = os.path.abspath(os.path.join(HERE, '..'))
FIG = os.path.join(NOTE, 'figures')

N = {}
for l in open(os.path.join(NOTE, 'logs', 'numbers.txt')):
    w = l.rstrip('\n').split(' ', 1)
    if len(w) == 2 and not w[0].startswith('CHECKS'):
        N[w[0]] = w[1].replace('{,}', ',')
def n(k): return N[k]
for k in ('ChkFill', 'ChkRepro', 'ChkTable', 'ChkInterp', 'ChkJobs', 'ChkEdge', 'ChkReproTwo', 'ChkEdgeTight', 'ChkKsSame'):
    assert N[k] == 'PASS', k
assert N['SampleStatus'].startswith('All 287'), 'note not in 287-segment mode'

d = Deck(date='Nov. 11, 2026', meeting='sPHENIX HF/Tracking meeting', figdir=FIG)
for sh in d.prs.slide_master.shapes:   # long meeting name: trim padding so "X.Dong/LBNL" stays on one line
    if sh.has_text_frame and 'Tracking' in sh.text_frame.text:
        r0 = sh.text_frame.paragraphs[0].runs[0]
        r0.text = re.sub(r'(meeting) +$', r'\1' + ' ' * 30, re.sub(r'(2026) +', r'\1' + ' ' * 40, r0.text))

K0 = 'K⁰S'; LA = 'Λ'; KAP = 'κ'

# 1 ------------------------------------------------------------------ title
d.title_slide('Per-Track Momentum Re-calibration (κ Table) Tested with K⁰S and Λ Mass Peaks',
              subtitle='sPHENIX run3pp p+p 200 GeV, KFParticle V0s, DDSim templates',
              lines=[f'{n("PerAllSeg")} segments, {n("NRuns")} runs, {n("PerAllEv")} events',
                     'κ(q, pT, η, φ) look-up table from runs 79509-79516 applied to every V0 daughter',
                     'Momentum scale and resolution refitted with and without κ',
                     'Tighter Λ selection for a cleaner Λ peak'])

# 2 ------------------------------------------------------------------ motivation
sl = d.slide('Motivation and Questions', layout='Title and Content', notes=
    'The K0S peak sits about 3 MeV below PDG and its width grows with pT much faster than in the full simulation. '
    'In the DDSim calibration this went into a global scale 0.991 and a curvature term 3.75 %/GeV, about 100 times the full '
    'simulation. A per-track correction table has been derived on runs 79509-79516. We apply it per daughter, recompute '
    'the masses and repeat the fit.')
d.bullets(sl, 0.4, 0.95, 9.2, 4.4, [
    f'run3pp {K0} peak ≈ 3 MeV/c² below PDG, width rises with pT much faster than in the full simulation',
    ('DDSim calibration absorbed it: s = 0.991, a = 3.75 %/(GeV/c), b = 1.40 % (current defaults)', {'lvl': 1}),
    'New input: per-track pT multiplier κ(q, pT, η, φ) derived on runs 79509-79516',
    (f'pT → κ·pT for each daughter; η, φ unchanged; masses recomputed', {'lvl': 1}),
    'Questions',
    ('Does κ remove the scale offset (s → 1)?', {'lvl': 1}),
    ('How much of the fitted resolution was a charge-, η-, φ-dependent scale distortion?', {'lvl': 1}),
    ('Does a table from runs 79509-79516 also work for later runs?', {'lvl': 1}),
    (f'Can a tighter {LA} selection pin down the {LA} peak (proton scale)?', {'lvl': 1}),
], size=17)
d.note_line(sl, 5.7, 'Same candidates, three masses: KFP fit / recomputed without κ / recomputed with κ')

# 3 ------------------------------------------------------------------ data
sl = d.slide('Data Set and Processing', notes=
    'Same 287 segments as the covariance-free production used for the DDSim calibration. Period A is the run range of the '
    'table, period B the other 24 runs where the table is applied as is. One condor job per segment reruns KFParticle and '
    'fills all three mass versions in one pass, so before and after never depend on two productions.')
d.bullets(sl, 0.4, 0.95, 9.2, 2.9, [
    'DST_TRKR_TRACKS, ana581_2026p005_v001, release ana.581 (streaming, MVTX misalignment, CDB)',
    'One condor job per segment: KFParticle K⁰S + Λ (no covariance cuts) → QA selection → mass vs pT',
    ('480714: KFP / no κ / with κ histograms, 287/287 OK, KFP histograms bin-identical to the earlier production', {'lvl': 1}),
    ('481293: + tight-Λ histograms + slim candidate trees (84 MB, cut studies now local), 287/287 OK, bin-identical', {'lvl': 1}),
    'Two periods: A = range of the κ table; B = all other runs (table applied as is)',
], size=15)
d.table(sl, 0.9, 4.05, 8.2, [
    ['period', 'runs', 'seg.', 'events', f'{K0} QA cand.', f'{LA} + anti-{LA} cand.'],
    ['A', f'79509-79516 ({n("NRunsA")})', n('PerASeg'), n('PerAEv'), n('PerANKs'), n('PerANLa')],
    ['B', f'{n("RunBFirst")}-{n("RunBLast")} ({n("NRunsB")})', n('PerBSeg'), n('PerBEv'), n('PerBNKs'), n('PerBNLa')],
], colw=[0.8, 2.1, 1.0, 1.2, 1.5, 1.6], size=13, rowh=0.36)
d.textbox(sl, 0.9, 5.25, 8.2, 0.4, ['QA candidates = all masses in the histogram range (signal + background)'], size=11, color=GREY)
d.note_line(sl, 5.75, 'All validation checks pass: fill macro, reproducibility, κ table, grid edges')

# 4 ------------------------------------------------------------------ selection
sl = d.slide('V0 Reconstruction and Selection', notes=
    'KFParticle pre-cuts are loose; all topological cuts are applied offline on stored values, identically for the three '
    'mass versions. No covariance, chi2 or significance cuts anywhere because DDSim cannot reproduce covariance correlations.')
d.table(sl, 0.3, 0.95, 9.4, [
    ['level', 'KFParticle (pre-cuts)', 'offline QA (data and DDSim templates)'],
    ['event', 'PV from SvtxVertexMap, crossing-matched', '|z_PV| < 10 cm'],
    ['track', 'pT > 0.1 GeV/c, ≥ 20 TPC hits, no χ²/DCA cut', '|η| < 1.1, |DCA_xy| > 0.02 cm (both daughters)'],
    ['pair', 'daughter DCA < 0.5 cm; vertex χ², volume off', 'DCA12 < 0.2 cm (SV-constrained → no effect)'],
    ['mother', 'DIRA > 0.8, mass window', 'DIRA > 0.995, L_xy > 0.5 cm'],
    ['species', 'K_S0 → π⁺π⁻, [Λ → pπ⁻]cc', 'Λ vs anti-Λ from the proton charge'],
    ['tight Λ', '—', '+ DIRA > 0.9995, |DCA_xy(π)| > 0.10, |DCA_xy(p)| > 0.03 cm'],
], colw=[1.0, 3.9, 4.5], size=12, rowh=0.42)
d.bullets(sl, 0.4, 4.15, 9.2, 2.2, [
    'No covariance-, χ²- or significance-based cut (DDSim cannot reproduce the covariance correlations)',
    'Cuts use the stored KFP values; they are not recomputed after κ',
    'All candidates kept (no best-candidate choice)',
], size=15)
d.note_line(sl, 5.75, 'Fits use the QA histograms; the tight Λ selection is a separate study (slides 13-14)')

# 5 ------------------------------------------------------------------ kappa table
sl = d.slide('The κ Look-up Table', notes=
    'Table convention from the file header: one over pT is multiplied by kappa_curv = 1 + eps + q delta, and kappa is its '
    'inverse. The charge-odd part dominates: negative tracks are scaled down, positive up, more at higher pT. Applied by '
    'linear interpolation in ln pT, eta, phi, clamped at the edges.')
d.bullets(sl, 0.4, 0.95, 5.4, 4.6, [
    f'{n("TabRows")} nodes: q = ±1 × {n("TabNpt")} pT ({n("TabPtMin")}-{n("TabPtMax")} GeV/c) × {n("TabNeta")} η × {n("TabNphi")} φ',
    'pT_corr = κ · pT,  κ = 1/(1 + ε + qδ)',
    ('ε: charge-even scale, qδ: charge-odd (sagitta-like) distortion', {'lvl': 1}),
    'Linear interpolation in (ln pT, η, φ), clamped at the edges',
    'Momentum vector scaled → η, φ unchanged; mass and V0 pT recomputed with PDG masses',
    f'Mean κ applied (QA, period A): {K0} π⁻ {n("PerAKapKsOne")}, π⁺ {n("PerAKapKsTwo")}; {LA} π {n("PerAKapLaOne")}, p {n("PerAKapLaTwo")}',
], size=15)
rows = [['pT (GeV/c)', 'κ (q = −1)', 'κ (q = +1)']]
for p, w in (('0.3', 'PtA'), ('1.0', 'PtB'), ('2.0', 'PtC')):
    rows.append([p, n(f'TabKMinus{w}'), n(f'TabKPlus{w}')])
d.table(sl, 6.1, 1.1, 3.5, rows, colw=[1.2, 1.15, 1.15], size=13, rowh=0.36)
d.textbox(sl, 6.1, 2.6, 3.5, 0.6, ['table nodes at η = 0, φ = 0'], size=11, color=GREY)
d.note_line(sl, 5.75, 'π⁺ and π⁻ are shifted in opposite directions → the K⁰S peak both moves and narrows')

# 6 ------------------------------------------------------------------ fit method
sl = d.slide('Momentum Scale / Resolution Fit', notes=
    'Same method as the DDSim main calibration. Templates are generated once without momentum smearing; for each grid '
    'point the daughter momenta are rescaled with fixed random numbers, so the chi2 surface is smooth. Per pT bin the data '
    'are fitted with N times the template plus a quadratic background, Poisson likelihood.')
d.bullets(sl, 0.4, 0.95, 9.2, 4.7, [
    'Response per daughter:  p_rec = s · p [1 + δ_p(pT) + σ G],   σ/pT = √((a pT)² + (b/β)²)',
    ('s scale, a curvature term, b multiple scattering; δ_p = full-sim proton bias (protons only)', {'lvl': 1}),
    'Templates: DDSim, 6M decays per species, no momentum smearing, same QA cuts; fixed random numbers',
    'Data: K⁰S in 7 V0-pT bins (0.2-4 GeV/c), Λ + anti-Λ in 5 bins (0.5-4 GeV/c, window 1.103-1.130 GeV/c²)',
    'Per bin: N × template + quadratic background, Poisson likelihood, Baker-Cousins χ² summed over bins',
    'Grid: a 0-8 %/(GeV/c), b 0.4-4 %, s 0.975-1.015 (widened: κ moves s to ≈ 1)',
    f'{K0}, {LA} and joint fits; Δχ²_Λ = {LA} χ² at the {K0} best point − {LA} minimum',
    'Periods A, B, all × KFP / no κ / with κ = 9 fits (+ 6 tight-Λ fits)',
], size=16)
d.note_line(sl, 5.75, '(a, b) are effective: they absorb every mass-width source DDSim does not model', color=RED)

# 7 ------------------------------------------------------------------ mass distributions
sl = d.slide('Mass Distributions Before and After κ', notes=
    'Same QA-selected candidates; blue recomputed without kappa, red with kappa, grey dashed the KFParticle mass. The K0S '
    'peak moves to the PDG mass and narrows in both periods. The Lambda peak narrows. The KFP mass has 20-25 percent more '
    'background entries in the same window, not understood.')
d.fig(sl, 'fig_massdist_int', 0.2, 0.85, w=6.6)
d.bullets(sl, 6.85, 1.0, 3.0, 5.3, [
    f'{K0} mean (MeV/c²)',
    (f'A {n("PerARawKsGM")} → {n("PerAKapKsGM")}', {'lvl': 1}),
    (f'B {n("PerBRawKsGM")} → {n("PerBKapKsGM")}', {'lvl': 1}),
    (f'PDG {n("PdgKs")}', {'lvl': 1}),
    f'{K0} σ (MeV/c²)',
    (f'A {n("PerARawKsGS")} → {n("PerAKapKsGS")}', {'lvl': 1}),
    (f'B {n("PerBRawKsGS")} → {n("PerBKapKsGS")}', {'lvl': 1}),
    f'{LA} σ (MeV/c²)',
    (f'A {n("PerARawLaGS")} → {n("PerAKapLaGS")}', {'lvl': 1}),
    (f'B {n("PerBRawLaGS")} → {n("PerBKapLaGS")}', {'lvl': 1}),
    ('Gaussian fits, pT-integrated', {'lvl': -1, 'color': GREY, 'size': 11}),
], size=14)
d.textbox(sl, 0.2, 5.75, 6.6, 0.4, ['blue: no κ, red: with κ, grey dashed: KFP mass; dotted: PDG'], size=11, color=GREY)
d.note_line(sl, 6.25, 'κ moves the K⁰S peak to within 0.4-0.6 MeV/c² of PDG and narrows both peaks')

# 8 ------------------------------------------------------------------ peaks vs pT
sl = d.slide('Peak Position and Width vs pT', notes=
    'Gaussian fits per V0-pT bin. Open: without kappa, full: with kappa; blue period A, red period B. The K0S width at high '
    'pT drops by about a third with kappa. Display only; the calibration uses the template fit.')
d.fig(sl, 'fig_peaks', 0.2, 0.85, w=6.3)
d.bullets(sl, 6.6, 1.0, 3.2, 5.0, [
    f'With κ the {K0} peak is within ≈ 1 MeV/c² of PDG in every pT bin (496.6-498.3)',
    f'{K0} width at V0 pT 1.6-2.2 GeV/c:',
    (f'A {n("PerARawKsWidHigh")} → {n("PerAKapKsWidHigh")} MeV/c²', {'lvl': 1}),
    (f'B {n("PerBRawKsWidHigh")} → {n("PerBKapKsWidHigh")} MeV/c²', {'lvl': 1}),
    f'{LA} peak stays ≈ 0.2-0.4 MeV/c² above PDG',
    ('open issue: proton scale', {'lvl': 1, 'color': RED}),
], size=14)
d.textbox(sl, 0.2, 5.5, 6.3, 0.5, ['open: no κ, full: with κ; blue A, red B; dashed PDG'], size=11, color=GREY)
d.note_line(sl, 6.2, 'Most of the high-pT width was a charge-dependent scale distortion')

# 9 ------------------------------------------------------------------ fit summary
sl = d.slide('Fit Results: Scale and Resolution', notes=
    'K0S fit; the joint fit gives the same points for the recomputed versions. Delta chi2 Lambda is the Lambda chi2 at the K0S '
    'point minus the Lambda minimum: it measures K0S-Lambda consistency.')
rows = [['period', 'mass', 'a (%/GeV)', 'b (%)', 's', 'χ²_K/389', 'Δχ²_Λ']]
for p, P in (('A', 'PerA'), ('B', 'PerB'), ('all', 'PerAll')):
    for v, V, lab in (('kfp', 'Kfp', 'KFP'), ('raw', 'Raw', 'no κ'), ('kappa', 'Kap', 'with κ')):
        rows.append([p if v == 'kfp' else '', lab, n(f'{P}{V}KsA'), n(f'{P}{V}KsB'), n(f'{P}{V}KsS'), n(f'{P}{V}KsChi'), n(f'{P}{V}DchiLa')])
d.table(sl, 0.6, 0.95, 8.8, rows, colw=[0.9, 1.3, 1.4, 1.1, 1.1, 1.5, 1.5], size=13, rowh=0.36)
d.bullets(sl, 0.4, 4.7, 9.2, 1.3, [
    f'Scale: s 0.991 → {n("PerAKapKsS")} (A), {n("PerBKapKsS")} (B); curvature term halved ({n("PerAARatio")} A, {n("PerBARatio")} B)',
    f'K⁰S-Λ tension in B removed by κ: Δχ²_Λ {n("PerBRawDchiLa")} → {n("PerBKapDchiLa")}',
], size=15)
d.note_line(sl, 6.1, f'Period B keeps a larger curvature term after κ ({n("PerBKapKsA")} vs {n("PerAKapKsA")} %/(GeV/c))', color=RED)

# 10 ----------------------------------------------------------------- resolution
sl = d.slide('Single-Track Resolution', notes=
    'Fitted sigma over pT at eta 0 for pions and protons, period A blue and B red, dashed without and solid with kappa; grey '
    'the full-simulation response. Bands are the grid nodes with Delta chi2 below 2.3; in period B only the best node passes.')
d.fig(sl, 'fig_res', 0.7, 0.85, w=8.6)
d.table(sl, 1.2, 5.05, 7.6, [
    ['π σ/pT (%)', '1 GeV/c', '2 GeV/c', '3 GeV/c'],
    ['A: no κ → κ', f'{n("PerARawResPiOne")} → {n("PerAKapResPiOne")}', f'{n("PerARawResPiTwo")} → {n("PerAKapResPiTwo")}', f'{n("PerARawResPiThree")} → {n("PerAKapResPiThree")}'],
    ['B: no κ → κ', f'{n("PerBRawResPiOne")} → {n("PerBKapResPiOne")}', f'{n("PerBRawResPiTwo")} → {n("PerBKapResPiTwo")}', f'{n("PerBRawResPiThree")} → {n("PerBKapResPiThree")}'],
], colw=[1.9, 1.9, 1.9, 1.9], size=13, rowh=0.33)
d.textbox(sl, 0.7, 4.66, 8.6, 0.35, ['blue A, red B; dashed no κ, solid with κ; grey: full simulation'], size=11, color=GREY)
d.note_line(sl, 6.3, 'With κ the pion resolution at 2 GeV/c roughly halves; still above the full simulation')

# 11 ----------------------------------------------------------------- chi2 maps
sl = d.slide('Δχ² Maps in (a, b)', notes=
    'K0S fit, minimised over s. Top period A, bottom B; left without, right with kappa. The best point moves to half the '
    'curvature term; the b term is stable. The minimum region is only a few grid nodes wide.')
d.fig(sl, 'fig_chi2', 0.2, 0.85, h=5.4)
d.bullets(sl, 7.25, 1.0, 2.6, 5.0, [
    'top: A, bottom: B',
    'left: no κ, right: with κ',
    'red star: best node',
    'a drops by 1.5 (A) and 2.0 (B) %/(GeV/c); b stays at 1.2-1.4 %',
    ('grid step (0.25 %/GeV in a) limits the interval resolution', {'lvl': 1}),
], size=14)
d.note_line(sl, 6.35, 'Well-separated minima: the κ effect is far larger than the statistical uncertainty')

# 12 ----------------------------------------------------------------- run dependence
sl = d.slide('Run Dependence', notes=
    'KFP-mass fits on run groups. The curvature term grows with the run period, and so does the number of K0S candidates '
    'per event, which points to luminosity or pile-up. Not tested here. The kappa table removes about half of a in both periods.')
rows = [['runs', 'a (%/(GeV/c)), KFP mass']]
rows += [['79509 seg. 0-30', n('GrpThirtyOne')], ['79509-79516 (A)', n('PerAKfpKsA')], ['80904-80915', n('GrpNine')],
         ['81556-81609', n('GrpEightOne')]]
d.table(sl, 0.6, 1.0, 5.0, rows, colw=[2.6, 2.4], size=14, rowh=0.38)
d.bullets(sl, 5.9, 1.0, 3.8, 3.2, [
    f'{K0} QA candidates per event:',
    (f'A {n("PerAKsEvMin")}-{n("PerAKsEvMax")}', {'lvl': 1}),
    (f'B {n("PerBKsEvMin")}-{n("PerBKsEvMax")}', {'lvl': 1}),
    'a grows with the same ordering',
], size=15)
d.bullets(sl, 0.4, 3.4, 9.2, 2.2, [
    'Hypothesis: luminosity-dependent distortion (e.g. TPC space charge) - not tested',
    f'κ (from period A) removes about half of a in both periods; B remains higher',
    'Candidate next step: κ table derived on period-B runs, or per-run fits',
], size=15)
d.note_line(sl, 5.75, 'Resolution is run dependent: quote DDSim parameters per period', color=RED)

# 13 ----------------------------------------------------------------- tight Lambda: choice
sl = d.slide('Tighter Λ Selection: Choice of Cuts', notes=
    'Scan on the 31 segments with full ntuples, side-band subtraction, only cuts DDSim models. DIRA is the strongest handle, '
    'decay length does not help because daughters need silicon. Bias test: DDSim selects before momentum smearing; if a cut '
    'picks well-measured candidates in data the fitted resolution would be biased. The width ratio tight over QA agrees.')
d.table(sl, 0.4, 0.95, 9.2, [
    ['31 segments, KFP mass, pT 0.5-4', 'S', 'B', 'S/B', 'S²/(S+B)', 'S kept'],
    ['QA selection', n('ScanBaseS'), n('ScanBaseB'), n('ScanBaseSB'), n('ScanBaseSeff'), '100 %'],
    ['+ DIRA > 0.9995', n('ScanDiraS'), n('ScanDiraB'), n('ScanDiraSB'), n('ScanDiraSeff'), n('ScanDiraKeep') + ' %'],
    ['+ L_xy > 2 cm (instead)', n('ScanLxyS'), n('ScanLxyB'), n('ScanLxySB'), n('ScanLxySeff'), n('ScanLxyKeep') + ' %'],
    ['tight: DIRA > 0.9995, |DCA(π)| > 0.10, |DCA(p)| > 0.03', n('ScanSelS'), n('ScanSelB'), n('ScanSelSB'), n('ScanSelSeff'), n('ScanSelKeep') + ' %'],
], colw=[4.2, 0.8, 0.8, 0.9, 1.3, 1.2], size=12, rowh=0.36)
d.bullets(sl, 0.4, 3.25, 9.2, 2.4, [
    f'{n("ScanNCuts")} combinations of DIRA, L_xy, π and p |DCA_xy|; only cuts DDSim can model',
    'Bias test (DDSim selects before smearing): peak width ratio tight / QA',
    (f'{K0}: data {n("BiasKsRatD")} ± {n("BiasKsRatDe")}, DDSim {n("BiasKsRatS")};  '
     f'{LA}: data {n("BiasLaRatD")} ± {n("BiasLaRatDe")}, DDSim {n("BiasLaRatS")}', {'lvl': 1}),
    (f'DDSim Λ templates keep {n("TmplLaTightFrac")} of their candidates', {'lvl': 1}),
], size=15)
d.note_line(sl, 5.75, f'Tight cut: S²/(S+B) × {n("ScanGain")}, no visible bias at the 5 % level (31 segments)')

# 14 ----------------------------------------------------------------- tight Lambda: results
sl = d.slide('Tighter Λ Selection: Results', notes=
    'Full sample, with kappa. Tight cuts: S over B up by four to five, peak-position error down by a quarter. Lambda-only '
    'scale agrees with K0S. Resolution: period A moves towards K0S, period B gets worse; the tight data widths at 1.2-2.4 GeV/c '
    'exceed the tight templates. So we use the tight selection for the scale only.')
d.fig(sl, 'fig_latight_mass', 0.2, 0.85, w=9.6)
rows = [['with κ', 'S/B', 'δ(peak) (MeV/c²)', 'S kept', 's_Λ', 'a_Λ (%/GeV)', 'Δχ²_Λ']]
for p, P in (('A', 'PerAKap'), ('B', 'PerBKap'), ('all', 'PerAllKap')):
    rows.append([f'{p}: QA → tight', f'{n(P+"LpQaSB")} → {n(P+"LpTiSB")}', f'{n(P+"LpQaMe")} → {n(P+"LpTiMe")}', n(P+'LpKeep') + ' %',
                 f'{n(P+"QaLaS")} → {n(P+"TiLaS")}', f'{n(P+"QaLaA")} → {n(P+"TiLaA")}', f'{n(P+"QaDchi")} → {n(P+"TiDchi")}'])
d.table(sl, 0.2, 4.65, 9.6, rows, colw=[1.45, 1.35, 1.55, 0.85, 1.55, 1.45, 1.4], size=11, rowh=0.3)
d.note_line(sl, 5.95, f'Used for the Λ scale (s_Λ = {n("PerAKapTiLaS")} A, {n("PerBKapTiLaS")} B), consistent with K⁰S', size=14)
d.note_line(sl, 6.35, 'Not for resolution: in B tight data widths exceed the tight templates (pointing modelling?)', color=RED, size=14)

# 15 ----------------------------------------------------------------- summary
sl = d.slide('Summary and Next Steps', layout='Title and Content', notes=
    'Summary. Decision needed on what DDSim should use.')
d.bullets(sl, 0.4, 0.95, 9.2, 4.0, [
    f'κ removes the momentum-scale offset: s 0.991 → {n("PerAllKapKsS")}; {K0} peak within 0.4-0.6 MeV/c² of PDG',
    f'Curvature term halved: a {n("PerAllRawKsA")} → {n("PerAllKapKsA")} %/(GeV/c) (all), b ≈ 1.4 % unchanged',
    f'Works outside its run range (period B) but B keeps a larger a ({n("PerBKapKsA")} vs {n("PerAKapKsA")})',
    f'{K0}/{LA} consistent with κ; tight Λ selection gives s_Λ = 0.999-1.000 (± 0.001)',
    ('Open', {'color': RED}),
    ('Λ peak 0.26 MeV/c² above PDG with κ: proton scale / δ_p double counting?', {'lvl': 1, 'color': RED}),
    ('Tight Λ widths in period B (DIRA vs DCA cuts, local on the slim trees)', {'lvl': 1, 'color': RED}),
    ('KFP vs recomputed mass: 20-25 % more background and different χ² in B', {'lvl': 1, 'color': RED}),
], size=16)
d.note_line(sl, 5.2, 'Decision: DDSim response = (a) uncorrected all-run defaults,', size=15)
d.note_line(sl, 5.55, '(b) κ-corrected per period (analyses apply κ), or (c) wait for a period-B table', size=15)
d.textbox(sl, 0.5, 6.1, 9.0, 0.4, ['Note: DataDrivenSim/kappa_note/build/main.pdf (28 pp)'], size=12, color=GREY)

# slide numbers
for i, s in enumerate(d.prs.slides, start=1):
    tb = s.shapes.add_textbox(Inches(9.15), Inches(7.12), Inches(0.7), Inches(0.32))
    pp = tb.text_frame.paragraphs[0]; pp.text = str(i); pp.alignment = PP_ALIGN.RIGHT
    for r in pp.runs: r.font.size = Pt(12); r.font.name = 'Times New Roman'; r.font.color.rgb = RGBColor(0, 0, 0)
d.save(os.path.join(HERE, 'kappa_momcal_summary.pptx'))
