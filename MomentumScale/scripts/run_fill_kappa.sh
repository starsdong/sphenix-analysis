#!/bin/bash
# run_fill_kappa.sh <mode>: V0 mass histograms from the run 79509 seg 0-30 loosecov KFP ntuples
#   mode 0 = KFP mass, 1 = recomputed from daughters, 2 = recomputed with the kappa pT correction
cd "$(dirname "$0")/.." || exit 1
K=data/loosecov/KFP_Ks_run79509_seg000-030_loosecov.root
L=data/loosecov/KFP_Lambda_run79509_seg000-030_loosecov.root
C=../MomentumScaleMap/run79509to79516_20261102/kappa_lookup.csv
mkdir -p data/kappa
root -l -b -q "macros/fillV0MassCorr.C+(\"$K\",\"$L\",\"$C\",\"data/kappa/V0mass_79509seg000-030_mode$1.root\",$1)" > data/kappa/fill_mode$1.log 2>&1
echo "mode $1 rc=$?"
