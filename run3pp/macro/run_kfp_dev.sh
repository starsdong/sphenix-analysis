#!/bin/bash
# DEV copy (uses macro/dev/Fun4All_run3pp_KFP.C, so production jobs are not affected by edits)
# usage: run_kfp_dev.sh <nEvents> <input DST or list> <outdir> [crossingStudy 0/1] [nSkip] [channels] [cutset loose|prod]
# Runs in a scratch dir, writes outputs to <outdir>, log to <outdir>/logs/.
nev=${1:-100}; input=$2; outdir=$3; xstudy=${4:-0}; nskip=${5:-0}; ch=${6:-Ks,Lambda,D0}; cutset=${7:-loose}
source /gpfs02/sphenix/user/dongxs/run3pp/setup_env.sh
# patched KFParticle_sPHENIX (single-track preselection for n>2-body decays), see src/coresoftware_kfp
export LD_LIBRARY_PATH=$TOP/install_kfpfix/lib:$LD_LIBRARY_PATH
[ "$xstudy" = "1" ] && xs=true || xs=false
mkdir -p $outdir/logs
base=$(basename $input .root)_$nskip
work=${_CONDOR_SCRATCH_DIR:-/tmp/dongxs_kfp_$$}; mkdir -p $work; cd $work
cp $TOP/macro/dev/Fun4All_run3pp_KFP.C $TOP/tools/slim_kfp.C .
log=$outdir/logs/$base.log
/usr/bin/time -v root.exe -b -q "Fun4All_run3pp_KFP.C($nev,\"$input\",\"$outdir\",$xs,$nskip,\"$ch\",\"$cutset\")" > $log 2>&1
rc=$?
if ! grep -q "^All done" $log; then echo "FAILED rc=$rc see $log"; exit 1; fi
echo "OK $log"
grep -m1 "^nEventsProcessed" $log
gzip -f $log   # logs are ~0.6 MB/job, mostly per-event lines
[ -z "$_CONDOR_SCRATCH_DIR" ] && cd / && rm -rf $work
exit 0
