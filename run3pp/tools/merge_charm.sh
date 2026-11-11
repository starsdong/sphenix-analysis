#!/bin/bash
# merge_charm.sh <production dir> : hadd the per-job KFP_<ch>_<run>_*.root files of a charm production per run and channel
# into <dir>/merged/KFP_<ch>_run<run>.root (nEventsProcessed is summed by hadd). Only jobs with a gzipped log (= OK) are used.
D=$(readlink -f $1)
source /gpfs02/sphenix/user/dongxs/run3pp/setup_env.sh > /dev/null 2>&1
cd $D && mkdir -p merged/lists
# outputs of successful jobs: first DST of the job list -> KFP_<ch>_<run>_<seg>_00000.root
ls -f logs | grep "\.txt_0\.log\.gz$" | sed "s/\.txt_0\.log\.gz$//" > merged/ok_jobs.txt
: > merged/ok_prefixes.txt
while read j; do head -1 lists/$j.txt | sed -E "s/.*-([0-9]{8})-([0-9]{5}).root/\1_\2/" >> merged/ok_prefixes.txt; done < merged/ok_jobs.txt
for ch in Dplus Ds Lc; do
  for r in $(cut -d_ -f1 merged/ok_prefixes.txt | sort -u); do
    grep "^${r}_" merged/ok_prefixes.txt | sed "s#^#$D/KFP_${ch}_#; s#\$#_00000.root#" > merged/lists/${ch}_${r}.txt
  done
done
ls merged/lists/*.txt | xargs -P 6 -I{} bash -c 'f={}; b=$(basename $f .txt); ch=${b%_*}; r=${b#*_}; hadd -f505 -k merged/KFP_${ch}_run${r}.root @$f > merged/lists/$b.hadd.log 2>&1 && echo "OK $b $(wc -l < $f)" || echo "FAILED $b"'
