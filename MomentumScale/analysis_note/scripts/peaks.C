// peaks.C -- Gaussian peak position / width of the Ks and Lambda+anti-Lambda mass per V0-pT bin for the two run periods
// and the three mass versions (kfp, raw, kappa). Prints GPEAK rows (machine-readable) and draws fig_peaks.
//   root -l -b -q 'scripts/peaks.C("<note dir>","<project dir>","<merged V0mass dir>")'
// Fits: Ks gaus + pol1 in 0.46-0.54, 2 MeV bins; Lambda gaus + pol2 in 1.103-1.130, 0.5 MeV bins. Mean limited to
// PDG +- 5 MeV, width 0.3-25 (Ks) or 0.3-10 (Lambda) MeV. A bin is drawn only if the fit is constrained (mean error < 2 MeV, width error < 60%).
#include <TFile.h>
#include <TH2F.h>
#include <TH1D.h>
#include <TF1.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TLine.h>
#include <TStyle.h>
#include <cstdio>
#include <cmath>
#include <vector>
#include <string>

namespace pk {
struct Res { double m, em, s, es, n; };
Res fit(TH2F* h, double plo, double phi, bool ks)
{
  TH1D* p = h->ProjectionY(Form("%s_%.1f_%.1f", h->GetName(), plo, phi), h->GetXaxis()->FindBin(plo + 1e-6), h->GetXaxis()->FindBin(phi - 1e-6));
  p->Rebin(2);
  const double m0 = ks ? 0.497611 : 1.115683, lo = ks ? 0.46 : 1.103, hi = ks ? 0.54 : 1.130;
  TF1 g("g", ks ? "gaus(0)+pol1(3)" : "gaus(0)+pol2(3)", lo, hi);
  g.SetParameters(p->GetMaximum(), m0, ks ? 0.005 : 0.0015, 0, 0, 0);
  g.SetParLimits(1, m0 - 0.005, m0 + 0.005); g.SetParLimits(2, 0.0003, ks ? 0.025 : 0.010);
  p->Fit(&g, "QNR");
  return {1e3 * g.GetParameter(1), 1e3 * g.GetParError(1), 1e3 * std::fabs(g.GetParameter(2)), 1e3 * g.GetParError(2),
          p->Integral(p->FindBin(lo + 1e-6), p->FindBin(hi - 1e-6))};
}
bool ok(const Res& r) { return r.em < 2.0 && r.es > 0 && r.es < 0.6 * r.s; }
}  // namespace pk

void peaks(const char* note, const char* proj, const char* dataDir)
{
  gStyle->SetOptStat(0); gStyle->SetOptTitle(0);
  gStyle->SetTitleSize(0.055, "XY"); gStyle->SetLabelSize(0.05, "XY"); gStyle->SetTitleOffset(1.45, "Y");
  const std::vector<double> ptK = {0.2, 0.5, 0.7, 0.9, 1.2, 1.6, 2.2, 4.0}, ptL = {0.5, 0.9, 1.2, 1.6, 2.4, 4.0};
  const char* per[2] = {"A", "B"};
  const char* ver[3] = {"kfp", "raw", "kappa"};
  TGraphErrors* g[2][2][3][2];   // channel, period, version, (mean, width)
  for (int c = 0; c < 2; ++c) {
    const bool ks = c == 0;
    const auto& pt = ks ? ptK : ptL;
    for (int p = 0; p < 2; ++p)
      for (int v = 0; v < 3; ++v) {
        TFile f(Form("%s/V0mass_%s_%s.root", dataDir, per[p], ver[v]));
        TH2F* h = (TH2F*)f.Get(ks ? "hKs_qa" : "hLa_qa")->Clone(Form("h%d%d%d", c, p, v)); h->SetDirectory(nullptr);
        if (!ks) h->Add((TH2F*)f.Get("hLaBar_qa"));
        for (int w = 0; w < 2; ++w) g[c][p][v][w] = new TGraphErrors();
        for (size_t j = 0; j + 1 <= pt.size(); ++j) {
          const bool all = j + 1 == pt.size();
          const double lo = all ? pt.front() : pt[j], hi = all ? pt.back() : pt[j + 1];
          pk::Res r = pk::fit(h, lo, hi, ks);
          printf("GPEAK %s %s %s %.1f %.1f %.3f %.3f %.3f %.3f %.0f %d\n", ks ? "Ks" : "Lambda", per[p], ver[v], lo, hi, r.m, r.em, r.s, r.es, r.n, int(pk::ok(r)));
          if (all || !pk::ok(r)) continue;
          const double x = 0.5 * (lo + hi) + (p == 0 ? -0.03 : 0.03), ex = 0.5 * (hi - lo);
          int n = g[c][p][v][0]->GetN();
          g[c][p][v][0]->SetPoint(n, x, r.m); g[c][p][v][0]->SetPointError(n, ex, r.em);
          g[c][p][v][1]->SetPoint(n, x, r.s); g[c][p][v][1]->SetPointError(n, ex, r.es);
        }
      }
  }
  TCanvas cv("cv", "", 1400, 1050);
  cv.Divide(2, 2);
  const char* yt[2] = {"peak position (MeV/c^{2})", "Gaussian #sigma (MeV/c^{2})"};
  for (int c = 0; c < 2; ++c)
    for (int w = 0; w < 2; ++w) {
      cv.cd(2 * c + w + 1);
      gPad->SetLeftMargin(0.17); gPad->SetBottomMargin(0.14); gPad->SetTopMargin(0.10); gPad->SetRightMargin(0.03);
      const double y0 = w ? 0 : (c ? 1112.5 : 490), y1 = w ? (c ? 6.0 : 22) : (c ? 1118.5 : 503);
      TH1* fr = gPad->DrawFrame(0, y0, 4.0, y1);
      fr->SetXTitle("V0 p_{T} (GeV/c)"); fr->SetYTitle(yt[w]);
      if (!w) { TLine* l = new TLine(0, c ? 1115.683 : 497.611, 4.0, c ? 1115.683 : 497.611); l->SetLineStyle(2); l->Draw(); }
      const double lx[4][4] = {{0.45, 0.17, 0.96, 0.36}, {0.20, 0.68, 0.70, 0.88}, {0.20, 0.16, 0.70, 0.35}, {0.20, 0.70, 0.70, 0.89}};
      const double* L = lx[2 * c + w];
      TLegend* lg = new TLegend(L[0], L[1], L[2], L[3]); lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.045);
      lg->SetNColumns(2);
      const int col[2] = {kBlue + 1, kRed + 1};
      for (int p = 0; p < 2; ++p)
        for (int v = 1; v < 3; ++v) {
          TGraphErrors* x = g[c][p][v][w];
          x->SetMarkerStyle(v == 1 ? 24 : 20); x->SetMarkerSize(1.3); x->SetMarkerColor(col[p]); x->SetLineColor(col[p]);
          x->Draw("P same");
          lg->AddEntry(x, Form("%s, %s", per[p], v == 1 ? "no #kappa" : "with #kappa"), "p");
        }
      lg->Draw();
      TLatex t; t.SetNDC(); t.SetTextSize(0.055);
      t.DrawLatex(0.18, 0.92, c ? "#Lambda+#bar{#Lambda}, QA cuts" : "K^{0}_{S}, QA cuts");
    }
  cv.SaveAs(Form("%s/figures/fig_peaks.pdf", note));
  cv.SaveAs(Form("%s/figures/fig_peaks.png", note));
}
