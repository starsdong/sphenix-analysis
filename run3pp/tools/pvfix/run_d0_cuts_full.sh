#!/bin/bash
# run_d0_cuts.sh : D0 search (pvfix_d0_v2.C) on the v2 skim with the pair-quality cuts, one output prefix per variant.
#   INTT  : both daughters have >= 1 INTT hit (crossing assignment from INTT timing)
#   COLL  : each daughter attached to the candidate's PV or |dz| < 1 mm to it (dz from the original-PV DCA)
cd /gpfs02/sphenix/user/dongxs/run3pp/output/pvfix2_10pct/analysis/full || exit 1
source ../../../../setup_env.sh > /dev/null 2>&1
run() { root -l -b -q "pvfix_d0_v2.C+(\"D0skim_pvfix2_10pct.root\",\"$1\",\"$2\",\"$3\")" > $1.log 2>&1; grep -E "rror" $1.log | head -3; echo "== $1"; cut -c1-260 $1.txt; }
run d0cut_none "1" "(no extra cut)"
run d0cut_intt "track_1_INTT_nHits>0 && track_2_INTT_nHits>0" "(INTT on both)"
run d0cut_coll "(track_1_inPV==1 || dz1<0.1) && (track_2_inPV==1 || dz2<0.1)" "(same collision)"
run d0cut_both "track_1_INTT_nHits>0 && track_2_INTT_nHits>0 && (track_1_inPV==1 || dz1<0.1) && (track_2_inPV==1 || dz2<0.1)" "(INTT + same collision)"
echo RUN_D0_CUTS_DONE
