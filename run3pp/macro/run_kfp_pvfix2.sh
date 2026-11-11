LD_LIBRARY_PATH now /sphenix/u/dongxs/install/lib:/cvmfs/sphenix.sdcc.bnl.gov/alma9.2-gcc-14.2.0/opt/sphenix/core/gcc/14.2.0-2f0a0/x86_64-el9/lib:/cvmfs/sphenix.sdcc.bnl.gov/alma9.2-gcc-14.2.0/opt/sphenix/core/gcc/14.2.0-2f0a0/x86_64-el9/lib64:/cvmfs/sphenix.sdcc.bnl.gov/alma9.2-gcc-14.2.0/opt/sphenix/core/binutils/2.40-acaab/x86_64-el9/lib:/cvmfs/sphenix.sdcc.bnl.gov/alma9.2-gcc-14.2.0/release/release_new/new.8/lib64:/cvmfs/sphenix.sdcc.bnl.gov/alma9.2-gcc-14.2.0/release/release_new/new.8/lib:/cvmfs/sphenix.sdcc.bnl.gov/alma9.2-gcc-14.2.0/opt/sphenix/utils/lib64:/cvmfs/sphenix.sdcc.bnl.gov/alma9.2-gcc-14.2.0/opt/sphenix/utils/lib:/cvmfs/sphenix.sdcc.bnl.gov/alma9.2-gcc-14.2.0/opt/sphenix/core/lib:/cvmfs/sphenix.sdcc.bnl.gov/alma9.2-gcc-14.2.0/opt/sphenix/core/root-6.32.06/lib:/cvmfs/sphenix.sdcc.bnl.gov/alma9.2-gcc-14.2.0/opt/sphenix/core/geant4.10.07.p04/lib64:/cvmfs/sphenix.sdcc.bnl.gov/alma9.2-gcc-14.2.0/opt/sphenix/core/lhapdf-5.9.1/lib:/usr/local/lib64:/usr/lib64
PATH now /sphenix/u/dongxs/install/bin:/cvmfs/sphenix.sdcc.bnl.gov/alma9.2-gcc-14.2.0/opt/sphenix/core/gcc/14.2.0-2f0a0/x86_64-el9/bin:/cvmfs/sphenix.sdcc.bnl.gov/alma9.2-gcc-14.2.0/opt/sphenix/core/binutils/2.40-acaab/x86_64-el9/bin:/cvmfs/sphenix.sdcc.bnl.gov/alma9.2-gcc-14.2.0/release/release_new/new.8/bin:/cvmfs/sphenix.sdcc.bnl.gov/alma9.2-gcc-14.2.0/opt/sphenix/utils/bin:/cvmfs/sphenix.sdcc.bnl.gov/alma9.2-gcc-14.2.0/opt/sphenix/core/bin:/cvmfs/sphenix.sdcc.bnl.gov/alma9.2-gcc-14.2.0/opt/sphenix/core/root-6.32.06/bin:/cvmfs/sphenix.sdcc.bnl.gov/alma9.2-gcc-14.2.0/opt/sphenix/core/geant4.10.07.p04/bin:/usr/local/bin:/usr/bin:/usr/local/sbin:/usr/sbin
#!/bin/bash
# D0 PV fixes: frozen macro in macro/pvfix2/, patched KFParticle in install_pvfix_v2/
# usage: run_kfp.sh <nEvents> <input DST or list> <outdir> [crossingStudy 0/1] [nSkip] [channels] [cutset loose|prod]
# Runs in a scratch dir, writes outputs to <outdir>, log to <outdir>/logs/.
nev=${1:-100}; input=$2; outdir=$3; xstudy=${4:-0}; nskip=${5:-0}; ch=${6:-Ks,Lambda,D0}; cutset=${7:-loose}
source /gpfs02/sphenix/user/dongxs/run3pp/setup_env.sh
# d0-pvfix: KFParticle_sPHENIX from src/coresoftware_kfp branch d0-pvfix (d0-pvfix head: cc077c6c5 3ba4ada8d c365374ff 59f846c13 4d032561d cb3844f6c 5c8c6eea5 ); post-selection via KFP_PVFIX_POSTCUT/KFP_PVFIX_LOOSE
export LD_LIBRARY_PATH=$TOP/install_pvfix_v2/lib:$LD_LIBRARY_PATH
export ROOT_INCLUDE_PATH=$TOP/install_pvfix_v2/include:$ROOT_INCLUDE_PATH
# production configuration of the d0-pvfix v2 run (overridable from the environment):
#   in-job D0 cuts loosened to track PV DCA >= 10 um, DIRA >= 0.5; production cuts (30 um, 0.9) applied at fill time
#   with the original PV OR the pair refit OR the KF refit w/o daughters; KF fit with track errors x7, chi2 cut 3.5
export KFP_PVFIX_LOOSE=${KFP_PVFIX_LOOSE:-0.001,0.5}
export KFP_PVFIX_POSTCUT=${KFP_PVFIX_POSTCUT:-0.003,0.9}
export KFP_PVKF_ERRSCALE=${KFP_PVKF_ERRSCALE:-7}
export KFP_PVKF_CHICUT=${KFP_PVKF_CHICUT:-3.5}
export KFP_PVFIX_PREFER_DIRA=${KFP_PVFIX_PREFER_DIRA:-0.9}  # keep the stock best combination whenever it passes the stock DIRA cut
export KFP_SLIM_DROP_ALLREFITS=${KFP_SLIM_DROP_ALLREFITS-1}  # in-job slimming drops *_PVrefitAll, *_PVkfAll (set to 0 to keep)
[ "$xstudy" = "1" ] && xs=true || xs=false
mkdir -p $outdir/logs
base=$(basename $input .root)_$nskip
work=${_CONDOR_SCRATCH_DIR:-/tmp/dongxs_kfp_$$}; mkdir -p $work; cd $work
cp $TOP/macro/pvfix2/Fun4All_run3pp_KFP.C $TOP/tools/slim_kfp.C .
log=$outdir/logs/$base.log
/usr/bin/time -v root.exe -b -q "Fun4All_run3pp_KFP.C($nev,\"$input\",\"$outdir\",$xs,$nskip,\"$ch\",\"$cutset\")" > $log 2>&1
rc=$?
if ! grep -q "^All done" $log; then echo "FAILED rc=$rc see $log"; exit 1; fi
echo "OK $log"
grep -m1 "^nEventsProcessed" $log
gzip -f $log   # logs are ~0.6 MB/job, mostly per-event lines
[ -z "$_CONDOR_SCRATCH_DIR" ] && cd / && rm -rf $work
exit 0
