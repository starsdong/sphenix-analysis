#!/bin/bash
# Wait for the 287-segment refits (k287_{B,all}_{kfp,raw,kappa}), then rebuild the kappa note in 287 mode.
#   bash kappa_note/scripts/watch_and_rebuild.sh [poll_s=120] [max_h=8]
# Exit 0: rebuilt; 1: a fit died without "All done"; 2: timeout; 3: build failed or still in 286 mode.
# Keeps the interim (286-seg) PDF and numbers as build/main_286seg.pdf, logs/numbers_286seg.txt and writes
# logs/numbers_diff_286_287.txt (every generated number that changed).
POLL=${1:-120}; MAXH=${2:-8}
NOTE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"; PROJ="$(cd "$NOTE/.." && pwd)"
M="$PROJ/output/momfit_kappa"; TAGS="B_kfp B_raw B_kappa all_kfp all_raw all_kappa"
t0=$(date +%s)
while :; do
  ndone=0; nrun=0
  for t in $TAGS; do
    grep -q "^All done" "$M/fit_k287_$t.log" 2>/dev/null && ndone=$((ndone+1))
    pgrep -f "fitMomRes.C.*\"k287_$t\"" >/dev/null && nrun=$((nrun+1))
  done
  echo "$(date '+%F %T') fits done $ndone/6, running $nrun"
  [ $ndone = 6 ] && [ $nrun = 0 ] && break
  if [ $((ndone+nrun)) -lt 6 ]; then
    echo "FAILED: a fit is neither running nor finished:"
    for t in $TAGS; do grep -q "^All done" "$M/fit_k287_$t.log" || pgrep -f "\"k287_$t\"" >/dev/null || { echo "  $t"; tail -5 "$M/fit_k287_$t.log"; }; done
    exit 1
  fi
  [ $(( $(date +%s) - t0 )) -gt $((MAXH*3600)) ] && { echo "TIMEOUT after ${MAXH} h"; exit 2; }
  sleep $POLL
done
cd "$NOTE"
[ -f build/main.pdf ] && cp -f build/main.pdf build/main_286seg.pdf
[ -f logs/numbers.txt ] && cp -f logs/numbers.txt logs/numbers_286seg.txt
bash build.sh > logs/rebuild_287.txt 2>&1; rc=$?
cat logs/rebuild_287.txt
grep -q "sample mode: 287" logs/rebuild_287.txt || { echo "BUILD NOT IN 287 MODE"; exit 3; }
[ $rc = 0 ] || { echo "BUILD FAILED rc=$rc"; exit 3; }
python3 - <<'E' > logs/numbers_diff_286_287.txt
def rd(f):
    d = {}
    for l in open(f):
        w = l.rstrip("\n").split(" ", 1)
        if len(w) == 2 and not w[0].startswith("CHECKS"): d[w[0]] = w[1]
    return d
a, b = rd("logs/numbers_286seg.txt"), rd("logs/numbers.txt")
for k in sorted(set(a) | set(b)):
    if a.get(k) != b.get(k): print(f"{k:22s} {a.get(k, '(new)')!s:>14s} -> {b.get(k, '(gone)')}")
E
echo "changed numbers: $(wc -l < logs/numbers_diff_286_287.txt)  (logs/numbers_diff_286_287.txt)"
grep -E '^PerB(Raw|Kap|Kfp)Ks[ABS] |^PerAll(Raw|Kap|Kfp)Ks[ABS] |DchiLa' logs/numbers.txt
printf -- "- %s: watch_and_rebuild.sh: 287-seg refits done, kappa_note rebuilt in 287 mode (%s); numbers diff in kappa_note/logs/numbers_diff_286_287.txt; hand-written text still to re-read.\n" \
  "$(date '+%F %H:%M %Z')" "$(grep -E '^pages' logs/rebuild_287.txt)" >> "$PROJ/session_log_2026-11-06.md"
echo "PDF: $NOTE/build/main.pdf (interim kept as build/main_286seg.pdf)"
