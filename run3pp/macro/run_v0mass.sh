#!/bin/bash
# usage: run_v0mass.sh <DST> <outdir>
# KFParticle Ks + Lambda without covariance cuts (macro/loosecov), then reduce to mass-vs-pT histograms
# (macro/loosecov/fill_v0_mass.C). Only V0mass_<base>.root and the gzipped logs are kept; the KFP ntuples stay in
# the job scratch dir and are deleted. Requires "All done" from both steps.
dst=$1; outdir=$2
TOPD=/gpfs02/sphenix/user/dongxs/run3pp
mkdir -p $outdir/logs
base=$(basename $dst .root)
work=${_CONDOR_SCRATCH_DIR:-/tmp/dongxs_v0mass_$$}; mkdir -p $work/kfp
bash $TOPD/macro/run_kfp_loosecov.sh 0 $dst $work/kfp 0 0 Ks,Lambda loose > $work/kfp_wrapper.out 2>&1
if ! grep -q "^OK" $work/kfp_wrapper.out; then echo "FAILED kfp"; cat $work/kfp_wrapper.out; cp $work/kfp/logs/* $outdir/logs/ 2>/dev/null; exit 1; fi
source $TOPD/setup_env.sh > /dev/null 2>&1
cp $TOPD/macro/loosecov/fill_v0_mass.C $work/
cd $work
root.exe -b -q "fill_v0_mass.C(\"$work/kfp\",\"$work/V0mass_$base.root\")" > $work/fill.log 2>&1
if ! grep -q "^All done" $work/fill.log; then echo "FAILED fill"; cat $work/fill.log; exit 1; fi
cp $work/V0mass_$base.root $outdir/V0mass_$base.root.tmp && mv $outdir/V0mass_$base.root.tmp $outdir/V0mass_$base.root
cp $work/kfp/logs/*.gz $outdir/logs/ 2>/dev/null
cat $work/kfp_wrapper.out $work/fill.log | grep -E "^OK|^nEventsProcessed|QA candidates"
rm -rf $work/kfp
[ -z "$_CONDOR_SCRATCH_DIR" ] && cd / && rm -rf $work
echo "OK $base"
exit 0
