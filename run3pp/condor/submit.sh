#!/bin/bash
# usage: submit.sh <tag> <listfile> <crossingStudy 0/1> [cutset loose|prod]
TOP=/gpfs02/sphenix/user/dongxs/run3pp
tag=$1; list=$(readlink -f $2); xs=${3:-0}; cs=${4:-loose}
out=$TOP/output/$tag; mkdir -p $out/logs /tmp/dongxs/condor_logs
sed -e "s#@OUTDIR@#$out#g" -e "s#@XSTUDY@#$xs#g" -e "s#@TAG@#$tag#g" -e "s#@CUTSET@#$cs#g" -e "s#@LIST@#$list#g" $TOP/condor/kfp.sub.template > $out/kfp.sub
if condor_q dongxs -af Cmd Args 2>/dev/null | grep -q "$out "; then echo "jobs for $out already queued"; exit 1; fi
condor_submit $out/kfp.sub
