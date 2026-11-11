#!/bin/bash
# run_d0_kc.sh <dir> <skim_kc> <prefix> : D0 search with kappa-corrected kinematics, no extra cut and INTT + same-collision
cd $1 || exit 1
source /gpfs02/sphenix/user/dongxs/run3pp/setup_env.sh > /dev/null 2>&1
cp /gpfs02/sphenix/user/dongxs/run3pp/output/pvfix2_10pct/analysis/pvfix_d0_v2_kc.C .
run() { root -l -b -q "pvfix_d0_v2_kc.C+(\"$2\",\"$3_$4\",\"$5\",\"$6\")" > $3_$4.log 2>&1; grep -E "rror" $3_$4.log | head -3; echo "== $3_$4"; cut -c1-260 $3_$4.txt; }
run $1 $2 $3 none "1" "(no extra cut)"
run $1 $2 $3 both "track_1_INTT_nHits>0 && track_2_INTT_nHits>0 && (track_1_inPV==1 || dz1<0.1) && (track_2_inPV==1 || dz2<0.1)" "(INTT + same collision)"
echo RUN_D0_KC_DONE
