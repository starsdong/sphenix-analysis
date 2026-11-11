#!/bin/bash
# merge_kappa287.sh: pull the condor 480714 per-segment outputs (run3pp/output/v0mass_kappa), merge per period and version,
# and print counts. Periods: A = runs 79509-79516 (range of the kappa table), B = all other runs (kappa applied as is).
# Output: data/kappa287/V0mass_{A,B,all}_{kfp,raw,kappa}.root
cd "$(dirname "$0")/.." || exit 1
R=/gpfs02/sphenix/user/dongxs/run3pp/output/v0mass_kappa
D=data/kappa287
mkdir -p $D/perjob
~/.hermes/scripts/sdcc/sx.sh "cd $R && tar czf ../v0mass_kappa_perjob.tgz V0mass_*.root && ls V0mass_*_kappa.root | wc -l" || exit 1
python3 ~/.hermes/scripts/sdcc/sxget.py /gpfs02/sphenix/user/dongxs/run3pp/output/v0mass_kappa_perjob.tgz $D || exit 1
tar xzf $D/v0mass_kappa_perjob.tgz -C $D/perjob
for v in kfp raw kappa; do
  A=$(ls $D/perjob/V0mass_*-0007950[9]-*_$v.root $D/perjob/V0mass_*-0007951[0-6]-*_$v.root)
  B=$(ls $D/perjob/V0mass_*_$v.root | grep -v -E -- '-00079509-|-0007951[0-6]-')
  echo "$v: period A $(echo $A | wc -w) files, period B $(echo $B | wc -w) files"
  hadd -f $D/V0mass_A_$v.root $A > /dev/null || exit 1
  hadd -f $D/V0mass_B_$v.root $B > /dev/null || exit 1
  hadd -f $D/V0mass_all_$v.root $D/perjob/V0mass_*_$v.root > /dev/null || exit 1
done
echo "merge done"
