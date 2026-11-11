// laTightPeaks.C -- Lambda+anti-Lambda Gaussian peak (gaus+pol2, 1.103-1.130, 0.5 MeV bins, likelihood) with QA vs tight cuts,
// per period and mass version, from data/kappa_v2 (condor 481293). Prints LPEAK rows; S, B = Gaussian integral and background
// within +-2 sigma.  root -l -b -q 'macros/laTightPeaks.C+'
#include <TFile.h>
#include <TH2.h>
#include <TH1D.h>
#include <TF1.h>
#include <cstdio>
#include <cmath>
void laTightPeaks()
{
  const char* per[3] = {"A", "B", "all"}; const char* ver[3] = {"kfp", "raw", "kappa"};
  const double bins[6][2] = {{0.5, 4.0}, {0.5, 0.9}, {0.9, 1.2}, {1.2, 1.6}, {1.6, 2.4}, {2.4, 4.0}};
  for (int p = 0; p < 3; ++p) for (int v = 0; v < 3; ++v) {
    TFile f(Form("data/kappa_v2/V0mass_%s_%s.root", per[p], ver[v]));
    for (int c = 0; c < 2; ++c) {
      TH2* h = (TH2*)f.Get(c ? "hLa_tight" : "hLa_qa")->Clone(Form("h%d%d%d", p, v, c)); h->Add((TH2*)f.Get(c ? "hLaBar_tight" : "hLaBar_qa"));
      for (auto& b : bins) {
        TH1D* x = h->ProjectionY(Form("x%d%d%d%.1f", p, v, c, b[0]), h->GetXaxis()->FindBin(b[0] + 1e-6), h->GetXaxis()->FindBin(b[1] - 1e-6));
        x->Rebin(2);
        TF1 g("g", "gaus(0)+pol2(3)", 1.103, 1.130);
        const double y0 = 0.5 * (x->GetBinContent(x->FindBin(1.104)) + x->GetBinContent(x->FindBin(1.129)));
        g.SetParameters(std::max(x->GetBinContent(x->FindBin(1.1157)) - y0, 1.), 1.1157, 0.0018, y0, 0, 0);
        g.SetParLimits(1, 1.1107, 1.1207); g.SetParLimits(2, 0.0003, 0.010);
        x->Fit(&g, "QNRL");
        const double m = g.GetParameter(1), s = std::fabs(g.GetParameter(2)), bw = x->GetBinWidth(1);
        const double S = g.GetParameter(0) * s * std::sqrt(2 * M_PI) / bw;
        TF1 bk("bk", "pol2", 1.103, 1.130); bk.SetParameters(g.GetParameter(3), g.GetParameter(4), g.GetParameter(5));
        const double B = bk.Integral(m - 2 * s, m + 2 * s) / bw;
        printf("LPEAK %s %s %s %.1f %.1f mean %.3f %.3f sigma %.3f %.3f S %.0f B2s %.0f SB %.2f\n", per[p], ver[v], c ? "tight" : "qa", b[0], b[1],
               1e3 * m, 1e3 * g.GetParError(1), 1e3 * s, 1e3 * g.GetParError(2), S, B, S / std::max(B, 1.));
      }
    }
  }
}
