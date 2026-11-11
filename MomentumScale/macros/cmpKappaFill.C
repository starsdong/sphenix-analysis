// cmpKappaFill.C -- bin-by-bin comparison of the cluster fill_v0_mass_kappa.C output with the local fillV0MassCorr.C output
// (same run 79509 seg 0-30 ntuples). Prints MAXDIFF per histogram and version; PASS if all are 0.
//   root -l -b -q 'macros/cmpKappaFill.C("<cluster prefix>","data/kappa/V0mass_79509seg000-030_mode")'
#include <TFile.h>
#include <TH2.h>
#include <cstdio>
#include <cmath>
void cmpKappaFill(const char* clusterPrefix, const char* localPrefix)
{
  const char* ver[3] = {"kfp", "raw", "kappa"};
  bool ok = true;
  for (int v = 0; v < 3; ++v) {
    TFile a(Form("%s_%s.root", clusterPrefix, ver[v])), b(Form("%s%d.root", localPrefix, v));
    for (const char* n : {"hKs_qa", "hKs_loose", "hLa_qa", "hLa_loose", "hLaBar_qa", "hLaBar_loose"}) {
      TH2* x = (TH2*)a.Get(n); TH2* y = (TH2*)b.Get(n);
      if (!x || !y) { printf("MISSING %s %s\n", ver[v], n); ok = false; continue; }
      double md = 0;
      for (int i = 0; i <= x->GetNbinsX() + 1; ++i)
        for (int j = 0; j <= x->GetNbinsY() + 1; ++j) md = std::max(md, std::fabs(x->GetBinContent(i, j) - y->GetBinContent(i, j)));
      // the local mode 0 file splits Lambda / anti-Lambda by the proton charge too, so all versions are comparable
      printf("MAXDIFF %-5s %-13s %g  (entries %.0f vs %.0f)\n", ver[v], n, md, x->GetEntries(), y->GetEntries());
      if (md > 0) ok = false;
    }
  }
  printf("CHECK %s\n", ok ? "PASS" : "FAIL");
}
