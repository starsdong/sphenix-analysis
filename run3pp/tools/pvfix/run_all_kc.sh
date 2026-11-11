#!/bin/bash
# run_all_kc.sh : after run_all_100.sh has finished (RUN_ALL_100_DONE in its log), apply the momentum-scale map to the
# full-sample D0 skim (kappa_apply.C -> D0skim_pvfix2_all_kc.root) and run the kappa-corrected D0 search
# (pvfix_d0_v2_kc.C) without extra cut and with INTT + same collision. Writes only into output/pvfix2_all/.
TOP=/gpfs02/sphenix/user/dongxs/run3pp
W=$TOP/output/pvfix2_all
cd $W || exit 1
while ! grep -q RUN_ALL_100_DONE run_all_100.log 2>/dev/null; do sleep 300; done
echo "$(date) run_all_100 done; applying the kappa map"
source $TOP/setup_env.sh > /dev/null 2>&1
root -l -b -q "$TOP/tools/pvfix/kappa_apply.C+(\"$TOP/calib/kappa_run79509to79516_20261102/kappa_lookup.csv\",\"D0skim_pvfix2_all.root\",\"D0skim_pvfix2_all_kc.root\")" 2>&1 | grep -E "kappa|D0_mass|rror"
sed 's/10%% sample/full sample/' $TOP/output/pvfix2_10pct/analysis/pvfix_d0_v2_kc.C > pvfix_d0_v2_kc.C
run() { root -l -b -q "pvfix_d0_v2_kc.C+(\"D0skim_pvfix2_all_kc.root\",\"$1\",\"$2\",\"$3\")" > $1.log 2>&1; grep -E "rror" $1.log | head -3; echo "== $1"; cut -c1-260 $1.txt; }
run d0allkc_none "1" "(no extra cut)"
run d0allkc_both "track_1_INTT_nHits>0 && track_2_INTT_nHits>0 && (track_1_inPV==1 || dz1<0.1) && (track_2_inPV==1 || dz2<0.1)" "(INTT + same collision)"
echo "$(date) RUN_ALL_KC_DONE"
