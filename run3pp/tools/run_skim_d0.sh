#!/bin/bash
# run_skim_d0.sh <file list> <output root file> : D0 second-stage skim (tools/skim_d0.C) in a condor job
source /gpfs02/sphenix/user/dongxs/run3pp/setup_env.sh > /dev/null 2>&1
work=${_CONDOR_SCRATCH_DIR:-/tmp/dongxs_skim_$$}
mkdir -p $work && cd $work
cp /gpfs02/sphenix/user/dongxs/run3pp/tools/skim_d0.C .
root -l -b -q "skim_d0.C+(\"$1\",\"$work/out.root\")" 2>&1 | grep -v -E "^Warning|^Info"
rc=${PIPESTATUS[0]}
if [ $rc -eq 0 ] && [ -s $work/out.root ]; then
  cp $work/out.root $2 && echo "OK $2"
else
  echo "FAILED rc=$rc"
  exit 1
fi
