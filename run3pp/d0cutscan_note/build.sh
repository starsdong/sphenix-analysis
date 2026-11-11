#!/bin/bash
# Build the D0 cut-scan note. Analysis outputs (logs/ tables/ figures/) come from scripts/d0cutscan.C, run with
#   bash build.sh --analysis      (re-runs all steps, ~15 min on the laptop, needs ../ntuples/d0skim and ../ntuples/d0tight_full)
set -e
ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$ROOT_DIR"
mkdir -p build logs tables figures
if [ "${1:-}" = "--analysis" ]; then
  source ~/root_v6/root_install/bin/thisroot.sh
  root-config --version > logs/root_version.txt
  for s in inputs scan grid trial ratio; do
    root -l -b -q -n "$ROOT_DIR/scripts/d0cutscan.C+(\"$s\")" > logs/run_$s.out 2>&1 || { echo "step $s FAILED"; tail -20 logs/run_$s.out; exit 1; }
  done
  for s in hump hump2 hump3; do   # origin of the pointing hump (Sec. 8b); hump.C figures are not used in the note
    root -l -b -q -n "$ROOT_DIR/scripts/$s.C+" > logs/run_$s.out 2>&1 || { echo "step $s FAILED"; tail -20 logs/run_$s.out; exit 1; }
  done
  rm -f figures/hump_theta.* figures/hump_predict.* figures/hump_slices.*
fi
python3 "$ROOT_DIR/scripts/make_numbers.py" "$ROOT_DIR"
for i in 1 2 3; do
  pdflatex -interaction=nonstopmode -halt-on-error -output-directory=build main.tex > build/pass$i.log 2>&1 || { echo "pdflatex pass $i FAILED"; tail -30 build/pass$i.log; exit 1; }
done
echo "undefined refs/citations: $(grep -c -E 'undefined' build/main.log || true)"
echo "missing files: $(grep -c -E 'not found' build/main.log || true)"
echo "overfull boxes: $(grep -c 'Overfull' build/main.log || true)"
used=$(grep -ho '\\includegraphics\[[^]]*\]{[^}]*}' sections/*.tex | sed 's/.*{\(.*\)}/\1/' | sort -u)
for f in $used; do [ -f figures/$f.pdf ] || [ -f figures/$f.png ] || echo "MISSING figure $f"; done
for f in $(ls figures/*.pdf | xargs -n1 basename | sed 's/\.pdf$//'); do echo "$used" | grep -qx "$f" || echo "unused figure $f"; done
echo "pages: $(pdfinfo build/main.pdf 2>/dev/null | awk '/^Pages/{print $2}')"
echo "PDF: $ROOT_DIR/build/main.pdf"
