// diagKappa.C -- machine-readable summary (DM, KAPPA, GPEAK rows) of at the kappa-corrected V0 mass files: recomputed-vs-KFP mass, mean kappa per daughter,
// mass peak position/width per V0-pT bin for mode 0/1/2 (Gaussian + pol1 in a window).
#include <TFile.h>
#include <TH1D.h>
#include <TH2F.h>
#include <TProfile.h>
#include <TF1.h>
#include <cstdio>
#include <cmath>
#include <vector>
#include <TString.h>
void diagKappa(const char* dir = "data/kappa")
{
  TString fs[3]; const char* f[3];
  for (int m = 0; m < 3; ++m) { fs[m] = Form("%s/V0mass_79509seg000-030_mode%d.root", dir, m); f[m] = fs[m].Data(); }
  TFile f1(f[1]), f2(f[2]);
  for (const char* h : {"hDm_Ks", "hDm_La"}) {
    TH1D* x = (TH1D*)f1.Get(h);
    printf("DM %s %.3f %.3f\n", h, 1e3 * x->GetMean(), 1e3 * x->GetRMS());   // recomputed - KFP mass (MeV): mean rms
  }
  for (const char* p : {"pKappa_Ks_d1", "pKappa_Ks_d2", "pKappa_La_d1", "pKappa_La_d2"}) {
    TProfile* x = (TProfile*)f2.Get(p);
    for (double pt : {0.15, 0.3, 0.5, 0.8, 1.2, 2.0, 3.0})
      printf("KAPPA %s %.2f %.4f %.0f\n", p, pt, x->GetBinContent(x->FindBin(pt)), x->GetBinEntries(x->FindBin(pt)));
  }
  struct Ch { const char* name; const char* h[2]; double m0, w, lo, hi; std::vector<double> pt; };
  Ch chs[2] = {{"Ks", {"hKs_qa", nullptr}, 0.4976, 0.005, 0.46, 0.54, {0.2, 0.5, 0.7, 0.9, 1.2, 1.6, 2.2, 4.0}},
               {"Lambda", {"hLa_qa", "hLaBar_qa"}, 1.1157, 0.002, 1.103, 1.130, {0.5, 0.9, 1.2, 1.6, 2.4, 4.0}}};
  for (auto& c : chs) {
    for (int m = 0; m < 3; ++m) {
      TFile F(f[m]);
      TH2F* h = (TH2F*)F.Get(c.h[0])->Clone("h"); if (c.h[1]) h->Add((TH2F*)F.Get(c.h[1]));
      for (size_t j = 0; j + 1 <= c.pt.size(); ++j) {
        const bool all = j + 1 == c.pt.size();
        const double lo = all ? c.pt.front() : c.pt[j], hi = all ? c.pt.back() : c.pt[j + 1];
        TH1D* p = h->ProjectionY(Form("p%d%zu", m, j), h->GetXaxis()->FindBin(lo + 1e-6), h->GetXaxis()->FindBin(hi - 1e-6));
        if (c.name[0] == 'K') p->Rebin(2); else p->Rebin(2);
        TF1 g("g", c.name[0] == 'K' ? "gaus(0)+pol1(3)" : "gaus(0)+pol2(3)", c.lo, c.hi);
        g.SetParameters(p->GetMaximum(), c.m0, c.w, 0, 0, 0); g.SetParLimits(2, 0.0003, 0.01); g.SetParLimits(1, c.m0 - 0.005, c.m0 + 0.005);
        p->Fit(&g, "QNR");
        printf("GPEAK %s %d %.1f %.1f %.3f %.3f %.3f %.3f %.0f\n", c.name, m, lo, hi, 1e3 * g.GetParameter(1), 1e3 * g.GetParError(1),
               1e3 * std::fabs(g.GetParameter(2)), 1e3 * g.GetParError(2), p->Integral());
      }
    }
  }
}
