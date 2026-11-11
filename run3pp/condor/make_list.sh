#!/bin/bash
# usage: make_list.sh <run> <first seg> <nfiles> > list   : DST paths of one run
D=/sphenix/lustre01/sphnxpro/production/run3pp/physics/ana581_2026p005_v001/DST_TRKR_TRACKS
run=$(printf %08d $1); lo=$(( $1/100*100 )); hi=$(( lo+100 ))
dir=$D/run_$(printf %08d $lo)_$(printf %08d $hi)
for s in $(seq $2 $(( $2+$3-1 ))); do f=$dir/DST_TRKR_TRACKS_run3pp_ana581_2026p005_v001-$run-$(printf %05d $s).root; [ -f $f ] && echo $f; done
