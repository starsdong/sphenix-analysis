#!/bin/bash
# watch_v2_and_fit.sh: wait for condor cluster 481293 (v0mass_kappa_v2), check all 287 jobs OK, merge
# (scripts/merge_kappa_v2.sh, incl. bin-by-bin reproducibility vs 480714), then run the tight-Lambda fits
# k287t_{A,B,all}_{raw,kappa} (fitMomRes laTight=1, s grid 0.975-1.015) and print the FIT lines.
#   bash scripts/watch_v2_and_fit.sh [poll_s=300] [max_h=16] [--skip-wait]
# Exit 0 done; 1 job failures/held; 2 timeout; 3 merge/repro failure; 4 fit failure.
POLL=${1:-300}; MAXH=${2:-16}; CL=481293
cd "$(dirname "$0")/.." || exit 1
R=/gpfs02/sphenix/user/dongxs/run3pp/output/v0mass_kappa_v2
t0=$(date +%s)
[ "$3" = "--skip-wait" ] && { echo "skipping the queue wait"; }
while [ "$3" != "--skip-wait" ]; do
  # an ssh failure (empty fields) is retried at the next poll
  st=$(~/.hermes/scripts/sdcc/sx.sh "condor_q $CL -af JobStatus 2>/dev/null | sort | uniq -c | tr '\n' ' '; echo; ls $R/V0mass_*_v0tree.root 2>/dev/null | wc -l; grep -l '^OK ' $R/logs/condor_*.out 2>/dev/null | wc -l; grep -l FAILED $R/logs/condor_*.out 2>/dev/null | wc -l" 2>/dev/null | tail -4)
  q=$(echo "$st" | sed -n 1p); ntree=$(echo "$st" | sed -n 2p | tr -d ' '); nok=$(echo "$st" | sed -n 3p | tr -d ' '); nfail=$(echo "$st" | sed -n 4p | tr -d ' ')
  echo "$(date '+%F %T') queue [${q:-empty}] trees $ntree OK $nok FAILED $nfail"
  if echo "$q" | grep -qE '(^| )5$|[0-9]+ 5( |$)'; then echo "HELD jobs in $CL: $q"; exit 1; fi
  [ "${nfail:-0}" -gt 0 ] && { echo "FAILED jobs: $nfail"; ~/.hermes/scripts/sdcc/sx.sh "grep -l FAILED $R/logs/condor_*.out | head"; exit 1; }
  if [ -z "$(echo "$q" | tr -d ' ')" ]; then
    [ "$nok" = 287 ] && [ "$ntree" = 287 ] && break
    echo "queue empty but OK=$nok trees=$ntree"; exit 1
  fi
  [ $(( $(date +%s) - t0 )) -gt $((MAXH*3600)) ] && { echo "TIMEOUT"; exit 2; }
  sleep $POLL
done
bash scripts/merge_kappa_v2.sh || exit 3
root -l -b -q -e 'gSystem->CompileMacro("macros/fitMomRes.C","k")' > /dev/null 2>&1
printf "%s\n" "k287t_A_raw data/kappa_v2/V0mass_A_raw.root" "k287t_A_kappa data/kappa_v2/V0mass_A_kappa.root" \
  "k287t_B_raw data/kappa_v2/V0mass_B_raw.root" "k287t_B_kappa data/kappa_v2/V0mass_B_kappa.root" \
  "k287t_all_raw data/kappa_v2/V0mass_all_raw.root" "k287t_all_kappa data/kappa_v2/V0mass_all_kappa.root" \
  | awk '{print $1, $2, "0.975 1.015 1"}' | xargs -P 6 -L 1 scripts/run_momfit_kappa.sh
ok=1
for t in A_raw A_kappa B_raw B_kappa all_raw all_kappa; do
  grep -q "^All done" output/momfit_kappa/fit_k287t_$t.log || { echo "fit k287t_$t did not finish"; ok=0; }
  grep -H -E "^TIGHT|^FIT |^CROSS" output/momfit_kappa/fit_k287t_$t.log
done
echo "- $(date '+%F %H:%M %Z'): condor $CL (v0mass_kappa_v2, tight Lambda + slim trees) done, merged to data/kappa_v2, tight fits k287t_* $( [ $ok = 1 ] && echo done || echo INCOMPLETE )." >> session_log_2026-11-06.md
[ $ok = 1 ] || exit 4
echo "ALL DONE"
