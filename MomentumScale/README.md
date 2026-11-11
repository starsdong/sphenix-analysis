# MomentumScale: track momentum scale / resolution from K0S and Lambda mass peaks (sPHENIX run3pp)

Calibration of the single-track momentum response sigma(pT)/pT = sqrt((a pT)^2 + (b/beta)^2) and scale s from the
K0S -> pi+ pi- and Lambda -> p pi- mass peaks in run3pp p+p 200 GeV data (287 segments, 32 runs, 285k events),
and test of a per-track pT re-calibration look-up table kappa(q, pT, eta, phi) (pT -> kappa * pT per daughter,
eta and phi unchanged). The fit uses templates from the data-driven fast simulation (DDSim).

Main results (K0S fit, without -> with kappa):

| period | a (%/(GeV/c)) | b (%) | s |
|---|---|---|---|
| A: runs 79509-79516 (range of the table) | 3.00 -> 1.50 | 1.60 -> 1.40 | 0.991 -> 0.999 |
| B: other 24 runs (table applied as is) | 4.00 -> 2.00 | 1.20 -> 1.40 | 0.991 -> 0.998 |
| all 287 segments | 3.50 -> 1.75 | 1.40 -> 1.40 | 0.991 -> 0.998 |

A tighter Lambda selection (DIRA > 0.9995, |DCA_xy(pi)| > 0.10 cm, |DCA_xy(p)| > 0.03 cm) raises the Lambda S/B
from 0.35 to 1.65 and gives a Lambda-only scale s = 0.999 (A) / 1.000 (B) with kappa; it is used for the scale only.
Statistical uncertainties only. Details: `analysis_note/kappa_momcal_note_v1.0.pdf`, summary slides in `slides/`.

## Layout

| directory | content |
|---|---|
| `cluster/` | SDCC condor production (run on sphnxuser, release ana.581) |
| `cluster/macro/loosecov/Fun4All_run3pp_KFP.C`, `cluster/macro/run_kfp_loosecov.sh` | KFParticle K0S / Lambda reconstruction without covariance cuts |
| `cluster/macro/kappa/fill_v0_mass_kappa.C` | candidate selection + mass-vs-pT histograms in three versions (KFP mass / recomputed / kappa-corrected) |
| `cluster/macro/kappa/fill_v0_mass_kappa_v2.C` | same + tight-Lambda histograms + slim V0 candidate tree |
| `cluster/macro/run_v0mass_kappa*.sh`, `cluster/condor/*.sub` | per-segment job wrapper (KFP in scratch, keep only histograms) and submit files |
| `cluster/macro/kappa/cmpKfpRepro.C`, `macros/cmpv2_cluster.C` | bin-by-bin reproducibility checks between productions |
| `macros/fitMomRes.C` | (a, b, s) grid fit: DDSim templates rescaled per grid point, per-pT-bin template + pol2 Poisson likelihood fit; option `laTight` |
| `macros/fillV0MassCorr.C`, `cmpKappaFill.C` | local implementation of the kappa fill and its validation against the cluster macro |
| `macros/diagKappa.C`, `plotKappaCompare.C`, `perRunV0.C` | diagnostics: kappa profiles, before/after comparison, per-run study |
| `macros/scanLambdaCuts.C`, `cutBiasTest.C`, `laTightPeaks.C` | tight-Lambda cut scan, cut-bias test, Lambda peaks QA vs tight |
| `scripts/` | fit launcher (`run_momfit_kappa.sh <tag> <file> [sLo sHi [laTight]]`), merge per run period, watcher |
| `analysis_note/` | LaTeX note (`main.tex`, `sections/`), figure/number scripts (`scripts/`), generated `tables/`, `figures/`, `logs/` (fit outputs the numbers are parsed from), built PDF |
| `slides/` | 15-slide summary (pptx + pdf) and its generator `make_slides.py` (numbers read from `analysis_note/logs/numbers.txt`) |

## Not included

- Data and ROOT files (per-segment histograms, slim trees, DDSim templates `Ks_nosmear_v2.root` / `Lambda_nosmear_v2.root`, fit outputs).
- The kappa look-up table itself (`MomentumScaleMap/run79509to79516_20261102/kappa_lookup.csv`, external input).
- The DDSim source (template generation); the templates were produced with `--no-mom-smear`, 6M decays per species.

`analysis_note/build.sh` and the scripts expect the original working-directory layout (DataDrivenSim project with
`output/momfit_kappa`, `data/kappa287`, `data/kappa_v2` next to the note); they document the procedure and are not
runnable from this repository alone.
