#!/bin/bash
# usage: run_v0mass_kappa_v2.sh <DST> <outdir>
# v2: fill_v0_mass_kappa_v2.C (v1 histograms + tight Lambda histograms + slim V0 candidate tree V0mass_<base>_v0tree.root)
# KFParticle Ks + Lambda without covariance cuts (macro/loosecov via run_kfp_loosecov.sh), then reduce to mass-vs-pT
# histograms in three versions (KFP mass / recomputed / kappa-corrected daughters) with macro/kappa/fill_v0_mass_kappa.C.
# Only V0mass_<base>_{kfp,raw,kappa,v0tree}.root and the gzipped logs are kept. Requires "All done" from both steps.
# Frozen inputs: macro/kappa/ (fill macro + kappa table) - do not edit while the production runs.
dst=$1; outdir=$2
TOPD=/gpfs02/sphenix/user/dongxs/run3pp
KCSV=$TOPD/macro/kappa/kappa_lookup_run79509to79516_20261102.csv
mkdir -p $outdir/logs
base=$(basename $dst .root)
work=${_CONDOR_SCRATCH_DIR:-/tmp/dongxs_v0kappa2_$$}; mkdir -p $work/kfp
bash $TOPD/macro/run_kfp_loosecov.sh 0 $dst $work/kfp 0 0 Ks,Lambda loose > $work/kfp_wrapper.out 2>&1
if ! grep -q "^OK" $work/kfp_wrapper.out; then echo "FAILED kfp"; cat $work/kfp_wrapper.out; cp $work/kfp/logs/* $outdir/logs/ 2>/dev/null; exit 1; fi
source $TOPD/setup_env.sh > /dev/null 2>&1
cp $TOPD/macro/kappa/fill_v0_mass_kappa_v2.C $KCSV $work/
cd $work
root.exe -b -q "fill_v0_mass_kappa_v2.C(\"$work/kfp\",\"$work/V0mass_$base\",\"$work/$(basename $KCSV)\")" > $work/fill.log 2>&1
if ! grep -q "^All done" $work/fill.log; then echo "FAILED fill"; cat $work/fill.log; exit 1; fi
for v in kfp raw kappa v0tree; do
  [ -s $work/V0mass_${base}_$v.root ] || { echo "FAILED missing $v"; exit 1; }
done
for v in kfp raw kappa v0tree; do
  cp $work/V0mass_${base}_$v.root $outdir/V0mass_${base}_$v.root.tmp && mv $outdir/V0mass_${base}_$v.root.tmp $outdir/V0mass_${base}_$v.root
done
cp $work/kfp/logs/*.gz $outdir/logs/ 2>/dev/null
cat $work/kfp_wrapper.out $work/fill.log | grep -E "^OK|^nEventsProcessed|QA candidates|kappa table|mean kappa|slim tree"
rm -rf $work/kfp
[ -z "$_CONDOR_SCRATCH_DIR" ] && cd / && rm -rf $work
echo "OK $base"
exit 0
