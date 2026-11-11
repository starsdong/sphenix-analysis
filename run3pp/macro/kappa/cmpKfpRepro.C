// cmpKfpRepro.C -- per segment: the kfp version of condor 480714 (output/v0mass_kappa/V0mass_<base>_kfp.root) must equal
// the condor 480608 file (output/v0mass_loosecov/V0mass_<base>.root) for hKs_qa, hKs_loose and the Lambda QA sum
// (480608 put Lambda+anti-Lambda in hLa_qa; 480714 splits them by the proton charge). Prints per-file DIFF lines and a summary.
//   root -l -b -q 'cmpKfpRepro.C("<new dir>","<old dir>")'
#include <TFile.h>
#include <TH2.h>
#include <TSystemDirectory.h>
#include <TList.h>
#include <TString.h>
#include <cstdio>
#include <cmath>
void cmpKfpRepro(const char* newDir, const char* oldDir)
{
  TSystemDirectory d("d", newDir);
  int n = 0, bad = 0, missing = 0;
  for (auto* o : *d.GetListOfFiles()) {
    TString f = o->GetName();
    if (!f.EndsWith("_kfp.root")) continue;
    TString base = f; base.ReplaceAll("_kfp.root", ""); base.ReplaceAll("V0mass_", "");
    TFile a(Form("%s/%s", newDir, f.Data())), b(Form("%s/V0mass_%s.root", oldDir, base.Data()));
    if (b.IsZombie()) { ++missing; printf("NOOLD %s\n", base.Data()); continue; }
    ++n;
    double md = 0;
    auto cmp = [&](TH2* x, TH2* y) {
      for (int i = 0; i <= x->GetNbinsX() + 1; ++i)
        for (int j = 0; j <= x->GetNbinsY() + 1; ++j) md = std::max(md, std::fabs(x->GetBinContent(i, j) - y->GetBinContent(i, j)));
    };
    cmp((TH2*)a.Get("hKs_qa"), (TH2*)b.Get("hKs_qa"));
    cmp((TH2*)a.Get("hKs_loose"), (TH2*)b.Get("hKs_loose"));
    TH2* la = (TH2*)a.Get("hLa_qa")->Clone("la"); la->Add((TH2*)a.Get("hLaBar_qa"));
    TH2* lb = (TH2*)b.Get("hLa_qa")->Clone("lb"); lb->Add((TH2*)b.Get("hLaBar_qa"));
    cmp(la, lb);
    if (md > 0) { ++bad; printf("DIFF %s maxdiff %g\n", base.Data(), md); }
  }
  printf("REPRO compared %d identical %d differ %d no-old %d\n", n, n - bad, bad, missing);
}
