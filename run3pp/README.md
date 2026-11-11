# run3pp: KFParticle charm reconstruction in sPHENIX Run-3 p+p

Code of the KFParticle study on the Run-3 p+p streaming track DSTs
(`ana581_2026p005_v001`, `DST_TRKR_TRACKS`, release `ana.581`): K0S / Lambda QA, D0 -> K pi with
primary-vertex (PV) bias fixes, momentum-scale (kappa) correction, and 3-prong charm
(D+ -> K pi pi, Ds+ -> phi pi, Lambda_c+ -> p K pi).

Code only. Data (ntuples, skims), figures, logs, result tables and built notes are not in this repository.

## Layout

| path | content |
|---|---|
| `setup_env.sh` | environment (release ana.581) |
| `macro/` | Fun4All macros and job wrappers. `Fun4All_run3pp_KFP.C` = stock production (Ks, Lambda, D0, 3-prong charm via `makeCharm3`); `pvfix/`, `pvfix2/` = frozen macros of the PV-fix productions; `run_kfp*.sh` = per-job wrappers |
| `condor/` | job-list making, submission templates, `check_prod.sh` (production status) |
| `src/EventInfoTree/` | small Fun4All module writing per-event info |
| `kfparticle_patches/` | `git format-patch` of branch `d0-pvfix` of coresoftware (7 commits, base `3ae2baa29`), the PV fixes in `KFParticle_sPHENIX` |
| `tools/` | ntuple QA, slimming (`slim_kfp.C`), skims, vertex-track audits, charm study tools |
| `tools/pvfix/` | PV-fix study: patch scripts `patch_fix*.py` (reproduce the commits), validation (`pvfix_sets.C`, `pvfix_validate.C`, `repro_v2.C`, `pv_shifts.C`, `pv_splithalf.C`, `trkstate_check.C`), pair studies, skims, kappa correction (`kappa_apply.C`, `kappa_v0check.C`), searches (`pvfix_d0_scan*.C`, `charm_search.C`), run scripts |
| `pvfix10pct/` | D0 search macros (`pvfix_d0_v2.C`, `_kc` = kappa-corrected, `_scan` = window scan) |
| `d0signal/` | earlier D0 signal / background studies |
| `analysis_note/`, `d0cutscan_note/`, `pvfix_note/` | LaTeX sources and number-generation scripts of the analysis notes (numbers are generated from logs by `scripts/make_numbers.py`; logs and figures not included) |
| `trackqa/` | track QA ntuple module, macro and note sources |

## PV fixes (branch `d0-pvfix`)

1. `primary_vertex_ID` = SvtxVertexMap key; `primary_vertex_crossing` added
2. PV refit without the candidate daughters (re-implemented PHSimpleVertexFinder pair average), `*_PVrefit`, `*_PVrefitAll`
3. KFParticle PV fit without daughters, `*_PVkf`, `*_PVkfAll`; 3b/3c: chi2 cut and track-error scale configurable
4. `pvfix_pass{Original,Refit,KF}` flags; production cuts applied with any PV at fill time (`KFP_PVFIX_POSTCUT`)
5. Best-combination choice restricted to stock-DIRA-passing combinations when one exists (stock candidates preserved)

Apply with `git am kfparticle_patches/*.patch` on coresoftware `3ae2baa29`, or run `tools/pvfix/patch_fix*.py`.
Run-time configuration: environment variables in `macro/run_kfp_pvfix2.sh`.

Paths in the scripts point to the SDCC areas used for the study (`/gpfs02/sphenix/user/dongxs/run3pp`).
