#!/bin/bash
# Build the kappa momentum-calibration note.   bash kappa_note/build.sh
# Copies fit logs from ../output/momfit_kappa, reruns the validation comparison, the figure macros and the numbers script,
# then pdflatex x3 and the standard checks.
set -e
NOTE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJ="$(cd "$NOTE/.." && pwd)"
KCSV="$PROJ/../MomentumScaleMap/run79509to79516_20261102/kappa_lookup.csv"
cd "$NOTE"
mkdir -p build logs tables figures
M="$PROJ/output/momfit_kappa"
# Sample mode: 287 when the refits with the last segment (run 81601 seg 05030, period B) are complete, else 286:
# period-B fits from output/momfit_kappa/k287_286seg, histograms from data/kappa287_286, no "all" fits.
MODE=287
for t in B_kfp B_raw B_kappa all_kfp all_raw all_kappa; do grep -q "^All done" "$M/fit_k287_$t.log" 2>/dev/null || MODE=286; done
echo "sample mode: $MODE segments"
DATA="$PROJ/data/kappa287"; [ $MODE = 286 ] && DATA="$PROJ/data/kappa287_286"
rm -rf fits; mkdir -p fits
TAGS="k287_A_kfp k287_A_raw k287_A_kappa k287_B_kfp k287_B_raw k287_B_kappa kfp31 g7951x g809xx g81xxx"
[ $MODE = 287 ] && TAGS="$TAGS k287_all_kfp k287_all_raw k287_all_kappa"
for t in $TAGS; do
  src="$M"; [ $MODE = 286 ] && case $t in k287_B_*) src="$M/k287_286seg";; esac
  grep -q "^All done" "$src/fit_$t.log" || { echo "fit $t not finished"; exit 1; }
  mkdir -p fits/$t; /bin/cp -f "$src/$t/"* fits/$t/
  /bin/cp -f "$src/$t/fitMomRes_$t.txt" logs/
  /bin/cp -f "$src/fit_$t.log" logs/
done
for p in A B; do for v in raw kappa; do for ch in Ks Lambda; do for e in png pdf; do
  /bin/cp -f "fits/k287_${p}_$v/fitMomRes_mass_${ch}_k287_${p}_$v.$e" figures/fig_mass_${ch}_${p}_$v.$e
done; done; done; done
# tight-Lambda study (condor 481293 = v2 production, data/kappa_v2; fits k287t_*; 31-seg scan and bias test)
for t in k287t_A_raw k287t_A_kappa k287t_B_raw k287t_B_kappa k287t_all_raw k287t_all_kappa; do
  grep -q "^All done" "$M/fit_$t.log" || { echo "tight fit $t not finished"; exit 1; }
  mkdir -p fits/$t; /bin/cp -f "$M/$t/"* fits/$t/
  /bin/cp -f "$M/$t/fitMomRes_$t.txt" logs/; /bin/cp -f "$M/fit_$t.log" logs/
done
for p in A B; do /bin/cp -f "fits/k287t_${p}_kappa/fitMomRes_mass_Lambda_k287t_${p}_kappa.png" figures/fig_mass_Lambda_${p}_tight.png; done
/bin/cp -f "$PROJ/output/scanLambdaCuts.log" logs/scanLambdaCuts.txt
/bin/cp -f "$PROJ/output/cutBiasTest.log" logs/cutBiasTest.txt
/bin/cp -f "$PROJ/data/kappa_v2/repro_v2.txt" logs/repro_v2.txt
(cd "$PROJ" && root -l -b -q 'macros/laTightPeaks.C+' 2>/dev/null | grep '^LPEAK') > logs/laTightPeaks.txt
root -l -b -q "scripts/latight.C(\"$NOTE\",\"$PROJ/data/kappa_v2\")" > /dev/null 2>&1
/bin/cp -f "$PROJ/analysis_note/logs/diagKappa.txt" logs/diagKappa.txt
/bin/cp -f "$PROJ/data/kappa/fill_mode2.log" logs/fill_validate.txt
/bin/cp -f "$PROJ/output/perRunV0.log" logs/perRunV0.txt
# cluster fill macro (run3pp/macro/kappa/fill_v0_mass_kappa.C) vs local fillV0MassCorr.C, run 79509 seg 0-30
(cd "$PROJ" && root -l -b -q 'macros/cmpKappaFill.C("data/kappa/cluster_validate/V0mass_79509seg000-030","data/kappa/V0mass_79509seg000-030_mode")' 2>/dev/null \
  | grep -E '^(MAXDIFF|CHECK|MISSING)') > logs/cmpKappaFill.txt
root -l -b -q "scripts/peaks.C(\"$NOTE\",\"$PROJ\",\"$DATA\")" 2>/dev/null | grep -E '^GPEAK' > logs/peaks.txt
root -l -b -q "scripts/massdist.C(\"$NOTE\",\"$DATA\")" 2>/dev/null | grep -E '^MDIST' > logs/massdist.txt
root -l -b -q "scripts/plot_res.C(\"$NOTE\",\"$PROJ\")" 2>/dev/null | grep -E '^(CURVE|BAND|MISSING)' > logs/res.txt
grep -q MISSING logs/res.txt && { echo "plot_res: missing fit"; exit 1; }
python3 scripts/make_numbers.py "$KCSV" $MODE > logs/numbers.txt
grep CHECKS logs/numbers.txt
for i in 1 2 3; do
  pdflatex -interaction=nonstopmode -halt-on-error -output-directory=build main.tex > build/pass$i.log 2>&1 || { echo "pdflatex pass $i FAILED"; tail -30 build/pass$i.log; exit 1; }
done
echo "undefined refs/citations: $(grep -c -E 'undefined' build/main.log || true)"
echo "missing files: $(grep -c -E 'not found' build/main.log || true)"
echo "overfull boxes: $(grep -c 'Overfull' build/main.log || true)"
echo "?? (unset macros) in text: $(pdftotext build/main.pdf - 2>/dev/null | grep -c '??' || true)"
used=$(grep -ho '\\includegraphics\[[^]]*\]{[^}]*}' sections/*.tex | sed 's/.*{\(.*\)}/\1/' | sort -u)
for f in $used; do [ -f figures/$f.pdf ] || [ -f figures/$f.png ] || echo "MISSING figure $f"; done
for f in $(ls figures/*.pdf 2>/dev/null | xargs -n1 basename | sed 's/\.pdf$//'); do echo "$used" | grep -qx "$f" || echo "unused figure $f"; done
echo "pages: $(pdfinfo build/main.pdf 2>/dev/null | awk '/^Pages/{print $2}')"
echo "PDF: $NOTE/build/main.pdf"
