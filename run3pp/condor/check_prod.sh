#!/bin/bash
# check_prod.sh <tag> <clusterId>  : status summary of a KFP condor production (run on an sphnxuser node)
# Scales to ~1e5 log files: no globbing / grepping over all logs (that hit "argument list too long" and was slow on GPFS).
#   done OK  = number of gzipped job logs (run_kfp.sh gzips the log only after "All done")
#   failed   = plain (not gzipped) job logs whose job is not in the queue any more
#   events   = sum of nEventsProcessed from condor_*.out, read incrementally (only newly finished jobs)
#   queued   = from the cluster ad (JobMaterializeNextProcId) because max_idle keeps jobs unmaterialized
# Appends a one-line record to output/<tag>/status_history.txt (used for rate / ETA).
TOP=/gpfs02/sphenix/user/dongxs/run3pp
tag=${1:-prod90pct_D0}; cl=${2:-480454}
out=$TOP/output/$tag
hist=$out/status_history.txt
st=$out/.check; mkdir -p $st
njobs=$(grep -c . $out/joblists.txt 2>/dev/null)
read idle run held < <(timeout 120 condor_q dongxs -constraint "ClusterId==$cl" -af JobStatus 2>/dev/null | awk '{c[$1]++} END{printf "%d %d %d\n", c[1]+0, c[2]+0, c[5]+0}')
next=$(timeout 60 condor_q $cl.-1 -af JobMaterializeNextProcId 2>/dev/null | head -1)
[[ "$next" =~ ^[0-9]+$ ]] || next=$njobs
waiting=$((njobs - next))
ls -f $out/logs > $st/ls.txt
grep '\.txt_0\.log\.gz$' $st/ls.txt | sed 's/\.txt_0\.log\.gz$//' | sort > $st/ok.txt
grep '\.txt_0\.log$' $st/ls.txt | sed 's/\.txt_0\.log$//' | sort > $st/plain.txt
ok=$(wc -l < $st/ok.txt)
# jobs still in the queue (their logs are plain while running)
timeout 120 condor_q dongxs -constraint "ClusterId==$cl" -af Args 2>/dev/null | awk '{n=split($2,a,"/"); sub(/\.txt$/,"",a[n]); print a[n]}' | sort > $st/queued.txt
comm -23 $st/plain.txt $st/queued.txt > $st/failed.txt
bad=$(wc -l < $st/failed.txt)
# incremental event count
touch $st/events.txt
cut -d' ' -f1 $st/events.txt | sort -u > $st/counted.txt   # resumes correctly after an interrupted run
comm -23 $st/ok.txt $st/counted.txt > $st/new.txt
while read j; do
  e=$(grep -m1 "^nEventsProcessed" $out/logs/condor_$j.out 2>/dev/null | awk '{print $2}')
  echo "$j ${e:-0}"
done < $st/new.txt >> $st/events.txt
cut -d' ' -f1 $st/events.txt | sort -u > $st/counted.txt
nevt=$(awk '!seen[$1]++{s+=$2} END{print s+0}' $st/events.txt)
nseg=$(sed "s#^#$out/lists/#; s#\$#.txt#" $st/ok.txt | xargs -r cat 2>/dev/null | grep -c .)
read nout size < <(find $out -maxdepth 1 -name "KFP_*.root" -printf "%s\n" 2>/dev/null | awk '{s+=$1; n++} END{printf "%d %.2f\n", n, s/1e9}')
now=$(date "+%Y-%m-%d %H:%M %Z"); ts=$(date +%s)
echo "[$tag cluster $cl] $now"
echo "  jobs: total $njobs | done OK $ok | FAILED $bad | running $run | idle $idle | held $held | not yet queued (max_idle) $waiting"
echo "  DST segments done: $nseg | events processed: $nevt | output files: $nout | size ${size} GB (running jobs' files are unslimmed)"
if [ -s $hist ]; then
  read t0 ok0 < <(tail -3 $hist | head -1 | awk '{for(i=1;i<=NF;i++){if($i~/^ts=/)t=substr($i,4); if($i~/^ok=/)o=substr($i,4)}; print t, o}')
  if [ -n "$t0" ] && [ $((ts - t0)) -gt 600 ] && [ $((ok - ok0)) -gt 0 ]; then
    awk -v d=$((ok - ok0)) -v dt=$((ts - t0)) -v rem=$((njobs - ok - bad)) \
      'BEGIN{r=d/dt*3600; printf "  rate %.0f jobs/h (recent) | ETA %.1f h at this rate\n", r, rem/r}'
  fi
fi
if [ "$held" -gt 0 ]; then
  echo "  held reasons:"; condor_q dongxs -constraint "ClusterId==$cl && JobStatus==5" -af ProcId HoldReason 2>/dev/null | cut -c1-150 | head -5 | sed 's/^/    /'
fi
if [ "$bad" -gt 0 ]; then
  echo "  failed jobs (first 10):"; head -10 $st/failed.txt | sed 's/^/    /'
fi
if [ $((idle + run + held)) -eq 0 ] && [ "$waiting" -eq 0 ] && [ "$njobs" -gt 0 ]; then
  echo "  >>> ALL JOBS LEFT THE QUEUE: $ok OK, $bad failed (of $njobs)."
fi
echo "$now ts=$ts idle=$idle run=$run held=$held ok=$ok bad=$bad seg=$nseg evt=$nevt size=$size" >> $hist
