#!/bin/bash
# merge_kappa_v2.sh: pull the condor 481293 per-segment outputs (run3pp/output/v0mass_kappa_v2: histograms incl. tight
# Lambda + slim V0 trees), check that the v1 histograms are bin-identical to production 480714, merge per period/version.
# Output: data/kappa_v2/V0mass_{A,B,all}_{kfp,raw,kappa}.root (hKs_qa, hLa_qa, hLaBar_qa, hLa_tight, hLaBar_tight, ...)
#         data/kappa_v2/v0tree/V0mass_<base>_v0tree.root (one per segment, for local cut studies)
cd "$(dirname "$0")/.." || exit 1
R=/gpfs02/sphenix/user/dongxs/run3pp/output
D=data/kappa_v2
mkdir -p $D/perjob $D/v0tree
# reproducibility vs 480714 on the cluster (v1 histograms of all three versions)
# (remote heredocs through sx.sh get mangled -> push the macro as a file)
S=$(ls -d ~/.hermes/skills/*/sphenix-sim-analysis/scripts | head -1)
(cd macros/cluster && bash $S/sxput.sh $R/v0mass_kappa_v2_validate cmpv2.C) || exit 1
~/.hermes/scripts/sdcc/sx.sh "cd $R && source ../setup_env.sh >/dev/null 2>&1 && root.exe -b -q v0mass_kappa_v2_validate/cmpv2.C 2>&1 | grep -E '^REPROV2|^DIFF' | tail -20" | tee $D/repro_v2.txt
grep -q "^REPROV2 segments 287 differing-histograms 0 missing-old 0" $D/repro_v2.txt || { echo "REPRO CHECK FAILED"; exit 1; }
~/.hermes/scripts/sdcc/sx.sh "cd $R/v0mass_kappa_v2 && tar cf ../v0mass_kappa_v2_perjob.tar V0mass_*.root && ls V0mass_*_v0tree.root | wc -l && du -sh ../v0mass_kappa_v2_perjob.tar" || exit 1
python3 ~/.hermes/scripts/sdcc/sxget.py $R/v0mass_kappa_v2_perjob.tar $D || exit 1
tar xf $D/v0mass_kappa_v2_perjob.tar -C $D/perjob && rm -f $D/v0mass_kappa_v2_perjob.tar
mv $D/perjob/*_v0tree.root $D/v0tree/
echo "trees: $(ls $D/v0tree | wc -l)"
for v in kfp raw kappa; do
  A=$(ls $D/perjob/V0mass_*-0007950[9]-*_$v.root $D/perjob/V0mass_*-0007951[0-6]-*_$v.root)
  B=$(ls $D/perjob/V0mass_*_$v.root | grep -v -E -- '-00079509-|-0007951[0-6]-')
  echo "$v: period A $(echo $A | wc -w) files, period B $(echo $B | wc -w) files"
  hadd -f $D/V0mass_A_$v.root $A > /dev/null || exit 1
  hadd -f $D/V0mass_B_$v.root $B > /dev/null || exit 1
  hadd -f $D/V0mass_all_$v.root $D/perjob/V0mass_*_$v.root > /dev/null || exit 1
done
echo "merge done"
