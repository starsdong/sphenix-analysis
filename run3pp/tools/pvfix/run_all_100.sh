#!/bin/bash
# run_all_100.sh : full-statistics (10% + 90%) analysis of the d0-pvfix v2 D0 production. Runs on sphnxuser01 under
# setsid/nohup; waits until condor cluster 480665 has left the queue, then:
#   1) reproduction check of the 90% against the stock 90% production (prod90pct_D0, same job lists): N(passOriginal)
#      per file must equal the stock candidate count
#   2) D0 skim of all finished OK jobs of both productions (skim_pvfix2.C)
#   3) D0 search (pvfix_d0_v2.C) without extra cuts, with INTT, with same-collision, with both
# Writes only into output/pvfix2_all/.
TOP=/gpfs02/sphenix/user/dongxs/run3pp
T=$TOP/tools/pvfix
W=$TOP/output/pvfix2_all
mkdir -p $W/done && cd $W || exit 1
while [ "$(condor_q 480665 -af ClusterId 2>/dev/null | wc -l)" != "0" ]; do sleep 600; done
echo "$(date) 480665 left the queue"
source $TOP/setup_env.sh > /dev/null 2>&1
python3 $T/done_files.py $TOP/output/pvfix2_10pct > done10.list
python3 $T/done_files.py $TOP/output/pvfix2_90pct > done90.list
echo "done files: 10% $(wc -l < done10.list) / $(grep -c . $TOP/output/pvfix2_10pct/joblists.txt), 90% $(wc -l < done90.list) / $(grep -c . $TOP/output/pvfix2_90pct/joblists.txt)"
for f in $(cat done10.list); do [ -e done/$f ] || ln -s $TOP/output/pvfix2_10pct/$f done/$f; done
for f in $(cat done90.list); do [ -e done/$f ] || ln -s $TOP/output/pvfix2_90pct/$f done/$f; done
echo "linked $(ls done | wc -l) (duplicate names across the two productions would show up as fewer links)"
root -l -b -q "$T/repro_v2.C+(\"$TOP/output/pvfix2_90pct\",\"$TOP/output/prod90pct_D0\",\"done90.list\")" 2>&1 | grep -E 'repro_v2|MISMATCH|rror' > repro90.txt
cat repro90.txt
root -l -b -q "$T/skim_pvfix2.C+(\"$W/done/KFP_D0_*.root\",\"D0skim_pvfix2_all.root\")" 2>&1 | grep -E 'skim_pvfix2|rror' > skim.txt
cat skim.txt
cp $TOP/output/pvfix2_10pct/analysis/pvfix_d0_v2.C .
sed -i 's/10%% sample/full sample/' pvfix_d0_v2.C
run() { root -l -b -q "pvfix_d0_v2.C+(\"D0skim_pvfix2_all.root\",\"$1\",\"$2\",\"$3\")" > $1.log 2>&1; grep -E "rror" $1.log | head -3; echo "== $1"; cut -c1-260 $1.txt; }
run d0all_none "1" "(no extra cut)"
run d0all_intt "track_1_INTT_nHits>0 && track_2_INTT_nHits>0" "(INTT on both)"
run d0all_coll "(track_1_inPV==1 || dz1<0.1) && (track_2_inPV==1 || dz2<0.1)" "(same collision)"
run d0all_both "track_1_INTT_nHits>0 && track_2_INTT_nHits>0 && (track_1_inPV==1 || dz1<0.1) && (track_2_inPV==1 || dz2<0.1)" "(INTT + same collision)"
echo "$(date) RUN_ALL_100_DONE"
