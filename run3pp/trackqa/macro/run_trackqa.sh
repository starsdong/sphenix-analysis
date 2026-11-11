#!/bin/bash
# usage: run_trackqa.sh <nEvents (<=0 all)> <input DST or list> <outdir> [trimDST 0/1] [minPt]
# Runs in a scratch dir; log to <outdir>/logs/. Requires "All done" in the log.
nev=${1:-0}; input=$2; outdir=$3; trim=${4:-0}; minpt=${5:-0}
source /gpfs02/sphenix/user/dongxs/run3pp/setup_env.sh
[ "$trim" = "1" ] && tr=true || tr=false
mkdir -p $outdir/logs
base=$(basename $input .root)
work=${_CONDOR_SCRATCH_DIR:-/tmp/dongxs_trackqa_$$}; mkdir -p $work; cd $work
cp $TOP/trackqa/macro/Fun4All_TrackQA.C .
log=$outdir/logs/$base.log
/usr/bin/time -v root.exe -b -q "Fun4All_TrackQA.C($nev,\"$input\",\"$outdir\",$tr,$minpt)" > $log 2>&1
rc=$?
if ! grep -q "^All done" $log; then echo "FAILED rc=$rc see $log"; exit 1; fi
echo "OK $log"
[ -z "$_CONDOR_SCRATCH_DIR" ] && cd / && rm -rf $work
exit 0
