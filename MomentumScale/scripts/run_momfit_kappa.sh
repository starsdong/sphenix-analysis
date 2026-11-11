#!/bin/bash
# run_momfit_kappa.sh <tag> <dataFile> [sLo sHi [laTight]]: (a, b, s) fit (macros/fitMomRes.C, main configuration: quadratic
# background, Lambda window 1.103-1.130, v2 no-smear templates) on a V0 mass file from fillV0MassCorr.C.
# Output: output/momfit_kappa/<tag>/, log output/momfit_kappa/fit_<tag>.log
cd "$(dirname "$0")/.." || exit 1
t=$1; d=$2; slo=${3:-0.984}; shi=${4:-1.000}; lt=${5:-0}
mkdir -p output/momfit_kappa
root -l -b -q "macros/fitMomRes.C+(\"$d\",\"output/momfit/Ks_nosmear_v2.root\",\"output/momfit/Lambda_nosmear_v2.root\",\"input/ddsim_input_run3pp_v2.root\",\"output/momfit_kappa/$t\",\"$t\",2000000,2,1.103,1.130,$slo,$shi,$lt)" > output/momfit_kappa/fit_$t.log 2>&1
echo "$t rc=$?"
