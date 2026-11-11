#!/bin/bash
# run_full_v2.sh [subdir] : analysis of the d0-pvfix v2 10% production on the jobs finished so far (or all, at the end).
# Uses only finished OK jobs (done_files.py) via symlinks in analysis/<subdir>/done/. Writes only into analysis/<subdir>/.
SUB=${1:-full}
TOP=/gpfs02/sphenix/user/dongxs/run3pp
O=$TOP/output/pvfix2_10pct
T=$TOP/tools/pvfix
mkdir -p $O/analysis/$SUB/done && cd $O/analysis/$SUB || exit 1
source $TOP/setup_env.sh > /dev/null 2>&1
python3 $T/done_files.py $O > done.list
for f in $(cat done.list); do [ -e done/$f ] || ln -s $O/$f done/$f; done
echo "done files: $(wc -l < done.list) of $(grep -c . $O/joblists.txt) job lists; linked $(ls done | wc -l)"
root -l -b -q "$T/repro_v2.C+(\"$O\",\"$TOP/output/pvfix10pct\",\"done.list\")" 2>&1 | grep -E 'repro_v2|MISMATCH|rror' > repro.txt
cat repro.txt
root -l -b -q "$T/skim_pvfix2.C+(\"$PWD/done/KFP_D0_*.root\",\"D0skim_pvfix2_10pct.root\")" 2>&1 | grep -E 'skim_pvfix2|rror' > skim.txt
cat skim.txt
cp $O/analysis/pvfix_d0_v2.C .
root -l -b -q "pvfix_d0_v2.C+(\"D0skim_pvfix2_10pct.root\",\"d0_pvfix2_10pct\")" > d0.log 2>&1
grep -E 'rror' d0.log | head -5
ls -la
echo RUN_FULL_V2_DONE
