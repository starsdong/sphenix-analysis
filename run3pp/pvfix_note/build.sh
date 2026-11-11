#!/bin/bash
# Build the d0-pvfix analysis note (10% preview). Inputs: logs/ and figures/ (copied from SDCC
# output/pvfix2_10pct/analysis/note_logs/, produced by tools/pvfix/note_logs.sh), tables from scripts/make_numbers.py.
set -e
ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$ROOT_DIR"
mkdir -p build logs tables figures
python3 "$ROOT_DIR/scripts/make_numbers.py" "$ROOT_DIR"
for i in 1 2 3; do
  pdflatex -interaction=nonstopmode -halt-on-error -output-directory=build main.tex > build/pass$i.log 2>&1 || { echo "pdflatex pass $i FAILED"; tail -30 build/pass$i.log; exit 1; }
done
echo "undefined refs/citations: $(grep -c -E 'undefined' build/main.log || true)"
echo "missing files: $(grep -c -E 'not found' build/main.log || true)"
echo "overfull boxes: $(grep -c 'Overfull' build/main.log || true)"
used=$(grep -ho '\\includegraphics\[[^]]*\]{[^}]*}' sections/*.tex | sed 's/.*{\(.*\)}/\1/' | sort -u)
for f in $used; do [ -f figures/$f.pdf ] || [ -f figures/$f.png ] || echo "MISSING figure $f"; done
echo "figures used: $(echo "$used" | grep -c .), available: $(ls figures/*.png | wc -l)"
echo "pages: $(pdfinfo build/main.pdf 2>/dev/null | awk '/^Pages/{print $2}')"
echo "PDF: $ROOT_DIR/build/main.pdf"
