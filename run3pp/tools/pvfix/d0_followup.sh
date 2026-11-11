#!/bin/bash
# d0_followup.sh <job> : full-sample D0 follow-up jobs (condor), output in output/pvfix2_all/followup/ (new files only)
#   0 kc, INTT+coll, window scan     1 kc, no cut, scan      2 uncorrected, INTT+coll, scan
#   3/4/5 kc, INTT+coll, runs 79509-79516 (kappa training) / 80880-80915 / 81556-81609
TOP=/gpfs02/sphenix/user/dongxs/run3pp
W=$TOP/output/pvfix2_all
O=$W/followup
mkdir -p $O && cd ${_CONDOR_SCRATCH_DIR:-$O} || exit 1
source $TOP/setup_env.sh > /dev/null 2>&1
cp $TOP/tools/pvfix/pvfix_d0_scan.C $TOP/tools/pvfix/pvfix_d0_scan_kc.C .
B="track_1_INTT_nHits>0 && track_2_INTT_nHits>0 && (track_1_inPV==1 || dz1<0.1) && (track_2_inPV==1 || dz2<0.1)"
KC=$W/D0skim_pvfix2_all_kc.root; UN=$W/D0skim_pvfix2_all.root
case $1 in
  0) m=pvfix_d0_scan_kc; in=$KC; p=fu_kc_both;  x="$B"; l="(full, INTT + same coll.)";;
  1) m=pvfix_d0_scan_kc; in=$KC; p=fu_kc_none;  x="1";  l="(full, no extra cut)";;
  2) m=pvfix_d0_scan;    in=$UN; p=fu_unc_both; x="$B"; l="(full, INTT + same coll.)";;
  3) m=pvfix_d0_scan_kc; in=$KC; p=fu_kc_both_r795; x="$B && runNumber>=79509 && runNumber<=79516"; l="(runs 79509-79516)";;
  4) m=pvfix_d0_scan_kc; in=$KC; p=fu_kc_both_r809; x="$B && runNumber>=80880 && runNumber<=80915"; l="(runs 80880-80915)";;
  5) m=pvfix_d0_scan_kc; in=$KC; p=fu_kc_both_r816; x="$B && runNumber>=81556 && runNumber<=81609"; l="(runs 81556-81609)";;
esac
root -l -b -q "$m.C+(\"$in\",\"$p\",\"$x\",\"$l\")" > $p.log 2>&1
cp $p.txt $p.log $p*_values.txt $p*.png $O/ 2>/dev/null
grep -q "^SCAN" $p.txt && echo "OK $p" || { echo "FAILED $p"; exit 1; }
