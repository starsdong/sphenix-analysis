#!/bin/bash
# run_charm.sh : charm_search.C on the full 10% stock-PV charm production (output/charm10pct), 3 channels in parallel,
# no extra cut and INTT on all daughters. Output: output/charm10pct/analysis/full/ (new files only).
TOP=/gpfs02/sphenix/user/dongxs/run3pp
A=$TOP/output/charm10pct/analysis/full
mkdir -p $A && cd $A || exit 1
source $TOP/setup_env.sh > /dev/null 2>&1
MAP=$TOP/calib/kappa_run79509to79516_20261102/kappa_lookup.csv
root -l -b -q -e ".L $TOP/tools/pvfix/charm_search.C+" > compile.log 2>&1
for ch in Dplus Ds Lc; do
  case $ch in Ds) it="phi_track_1_INTT_nHits>0 && phi_track_2_INTT_nHits>0 && track_3_INTT_nHits>0";;
              *) it="track_1_INTT_nHits>0 && track_2_INTT_nHits>0 && track_3_INTT_nHits>0";; esac
  ( root -l -b -q "$TOP/tools/pvfix/charm_search.C+(\"$ch\",\"$TOP/output/charm10pct/KFP_${ch}_*.root\",\"$MAP\",\"c_${ch}_none\",\"1\",\"(no extra cut)\")" > c_${ch}_none.log 2>&1
    root -l -b -q "$TOP/tools/pvfix/charm_search.C+(\"$ch\",\"$TOP/output/charm10pct/KFP_${ch}_*.root\",\"$MAP\",\"c_${ch}_intt\",\"$it\",\"(INTT)\")" > c_${ch}_intt.log 2>&1
    echo "$(date) $ch done" ) &
done
wait
echo "$(date) RUN_CHARM_DONE"
