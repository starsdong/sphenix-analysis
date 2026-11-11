#!/bin/bash
# note_logs.sh : collect / regenerate every input of the pvfix analysis note (10% preview) into
# output/pvfix2_10pct/analysis/note_logs/ (new files only). Run on sphnxuser01.
TOP=/gpfs02/sphenix/user/dongxs/run3pp
T=$TOP/tools/pvfix
O=$TOP/output
L=$O/pvfix2_10pct/analysis/note_logs
mkdir -p $L && cd $L || exit 1
source $TOP/setup_env.sh > /dev/null 2>&1
echo "== software"
{ echo "release: $(grep -m1 -o 'ana\.[0-9]*' $TOP/setup_env.sh)"; echo "OFFLINE_MAIN $OFFLINE_MAIN"; root-config --version | sed 's/^/root /';
  cd $TOP/src/coresoftware_kfp && echo "branch $(git rev-parse --abbrev-ref HEAD)"; git log --format='commit %h %s' -7 d0-pvfix; cd $L; } > software.txt 2>&1
md5sum $TOP/install_pvfix_v2/lib/libkfparticle_sphenix.so $TOP/calib/kappa_run79509to79516_20261102/kappa_lookup.csv >> software.txt
cat $TOP/macro/run_kfp_pvfix2.sh > run_kfp_pvfix2.sh.txt
echo "== jobs"
$TOP/condor/check_prod.sh pvfix2_10pct 480634 > jobs_v2_10pct.txt 2>&1
$TOP/condor/check_prod.sh pvfix10pct 480622 > jobs_v1_10pct.txt 2>&1
cp $O/pvfix2_10pct/kfp.sub kfp_v2_10pct.sub.txt
echo "== smoke"
S=$O/pvfix_smoke
F1=KFP_D0_00079513_00010_00000.root; F2=$(basename $(ls $S/stock_81609/*.root))
{ echo "SMOKE run79513"; root -l -b -q "$T/pvfix_sets.C+(\"$S/stock/$F1\",\"$S/v2fin/$F1\",\"$S/v2prod/$F1\")" 2>&1 | grep -E "^A:|^B:|^A not";
  echo "SMOKE run81609"; root -l -b -q "$T/pvfix_sets.C+(\"$S/stock_81609/$F2\",\"$S/v2fin_81609/$F2\",\"$S/v2prod_81609/$F2\")" 2>&1 | grep -E "^A:|^B:|^A not";
  echo "VALIDATE stock closure v2fin run79513"; root -l -b -q "$T/pvfix_validate.C+(\"$S/stock/$F1\",\"$S/closure/$F1\",\"$S/v2fin/$F1\")" 2>&1 | grep -E "^candidates|differs|common branches|^fix 1";
  for c in stock v2fin stock_81609 v2fin_81609 v2slim v2slim_81609; do
    f=$(ls $S/$c/*.root); echo "CPU $c user $(zcat $S/$c/logs/*.gz | grep 'User time' | awk '{print $4}') size $(stat -c %s $f) entries $(root -l -b -q -e "TFile f(\"$f\"); printf(\"%lld %d\n\", ((TTree*)f.Get(\"DecayTree\"))->GetEntries(), ((TTree*)f.Get(\"DecayTree\"))->GetListOfBranches()->GetEntries());" 2>/dev/null | tail -1)"; done; } > smoke.txt 2>&1
echo "== DST track states / split half"
root -l -b -q "$T/trkstate_check.C+(\"$(head -1 $O/vtxaudit/test.list)\",200)" > trkstate.txt 2>&1
cp $O/vtxaudit/split.log split_s1.txt; cp $O/vtxaudit/split_s7.log split_s7.txt; cp $O/vtxaudit/split_s10.log split_s10.txt
echo "== PV versions on 400 v2 files"
ls $O/pvfix2_10pct/analysis/full/done/*.root | awk 'NR%13==0' | head -400 > shift_files.list
mkdir -p shiftdir; for f in $(cat shift_files.list); do ln -sf $f shiftdir/; done
root -l -b -q "$T/pv_shifts.C+(\"$L/shiftdir/KFP_D0_*.root\")" > pv_shifts.txt 2>&1
echo "== pair crossing / pair cuts"
for r in 00079513 00081609; do root -l -b -q "$T/pair_crossing.C+(\"$O/pvfix2_10pct/KFP_D0_${r}_00000_00000.root\")" > pair_crossing_$r.txt 2>&1; done
root -l -b -q "$T/pair_cuts_study.C+(\"$O/pvfix2_10pct/analysis/full/D0skim_pvfix2_10pct.root\")" > pair_cuts.txt 2>&1
echo "== kappa"
cp $O/kappa_check/v0check.txt kappa_v0check.txt; cp $O/kappa_check/v0check_K0S.png $O/kappa_check/v0check_Lambda.png .
root -l -b -q "$T/kc_stats.C+(\"$O/pvfix2_10pct/analysis/full/D0skim_pvfix2_10pct_kc.root\")" > kc_stats.txt 2>&1
echo "== analysis outputs"
A=$O/pvfix2_10pct/analysis/full
cp $A/repro.txt $A/skim.txt .
for v in d0cut_none d0cut_intt d0cut_coll d0cut_both d0kc_none d0kc_both; do cp $A/$v.txt $A/${v}_values.txt .; for s in HiP HiT AllP AllT; do cp $A/${v}_$s.png .; done; done
echo NOTE_LOGS_DONE
