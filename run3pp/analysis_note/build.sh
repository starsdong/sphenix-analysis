#!/bin/bash
# Build the analysis note. Numbers/figures come from scripts/make_numbers.C (run on the cluster,
# outputs copied to logs/ tables/ figures/); this script regenerates job stats and the PDF.
set -e
ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$ROOT_DIR"
mkdir -p build
python3 "$ROOT_DIR/scripts/job_stats.py" "$ROOT_DIR" > /dev/null
for i in 1 2 3; do
  pdflatex -interaction=nonstopmode -halt-on-error -output-directory=build main.tex > build/pass$i.log 2>&1 || { echo "pdflatex pass $i FAILED"; tail -30 build/pass$i.log; exit 1; }
done
echo "undefined refs/citations: $(grep -c -E 'undefined' build/main.log || true)"
echo "missing files: $(grep -c -E 'not found|File .* not found' build/main.log || true)"
echo "overfull boxes: $(grep -c 'Overfull' build/main.log || true)"
# figures referenced vs present
used=$(grep -ho '\\includegraphics\[[^]]*\]{[^}]*}' sections/*.tex | sed 's/.*{\(.*\)}/\1/' | sort -u)
for f in $used; do [ -f figures/$f.pdf ] || [ -f figures/$f.png ] || echo "MISSING figure $f"; done
for f in $(ls figures/*.pdf | xargs -n1 basename | sed 's/\.pdf$//'); do echo "$used" | grep -qx "$f" || echo "unused figure $f"; done
echo "pages: $(pdfinfo build/main.pdf 2>/dev/null | awk '/^Pages/{print $2}')"
echo "PDF: $ROOT_DIR/build/main.pdf"
