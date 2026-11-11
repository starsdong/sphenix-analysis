// plot_res.C -- fitted single-track resolution sigma(pT)/pT at eta = 0 for pions and protons, run periods A and B,
// without and with kappa, plus the full-simulation input; and the Delta chi2 (a, b) maps (min over s) of the Ks fits.
//   root -l -b -q 'scripts/plot_res.C("<note dir>","<project dir>")'
// Best points: "FIT Ks" lines (Ks fit, = summary table); bands from chi2_ab_Ks; of logs/fitMomRes_k287_<P>_<v>.txt; maps: <note>/fits/k287_<P>_<v>/fitMomRes_*.root (staged by build.sh).
// Prints CURVE rows (sigma/pT at 1, 2, 3 GeV/c) and BAND rows (Delta chi2 < 2.3 envelope at 2 GeV/c, pions).
#include <TFile.h>
#include <TH2D.h>
#include <TH1.h>
#include <TF1.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TStyle.h>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>
#include <cstdio>

namespace pr {
std::vector<double> best(const std::string& log, const std::string& which)
{
  std::ifstream in(log);
  std::string line;
  while (std::getline(in, line)) {
    std::istringstream is(line);
    std::string w, name, x;
    is >> w >> name;
    if (w != "FIT" || name != which) continue;
    double a, b, s, d;
    is >> x >> a >> d >> d >> x >> b >> d >> d >> x >> s;
    return {a, b, s};
  }
  return {};
}
double res(double a, double b, double pt, double m)
{
  const double beta = pt / std::sqrt(pt * pt + m * m);
  return std::sqrt(a * a * pt * pt + b * b / (beta * beta));
}
}  // namespace pr

void plot_res(const char* note, const char* proj)
{
  gStyle->SetOptStat(0); gStyle->SetOptTitle(0);
  struct F { const char* tag; const char* label; int col; int lsty; };
  std::vector<F> fits = {{"k287_A_raw", "A, no #kappa", kBlue + 1, 2}, {"k287_A_kappa", "A, with #kappa", kBlue + 1, 1},
                         {"k287_B_raw", "B, no #kappa", kRed + 1, 2}, {"k287_B_kappa", "B, with #kappa", kRed + 1, 1}};
  const char* sp[2] = {"pi", "p"};
  const double m[2] = {0.13957, 0.93827};
  const char* lab[2] = {"#pi^{#pm}", "p, #bar{p}"};
  TFile fi(Form("%s/input/ddsim_input_run3pp_v2.root", proj));
  const int n = 120;
  TCanvas c("c", "", 1300, 600);
  c.Divide(2, 1);
  for (int i = 0; i < 2; ++i) {
    c.cd(i + 1);
    gPad->SetLeftMargin(0.16); gPad->SetBottomMargin(0.13); gPad->SetTopMargin(0.04); gPad->SetRightMargin(0.04);
    TH1* fr = gPad->DrawFrame(0.1, 0, 3.0, 0.26);
    fr->GetXaxis()->SetTitle("p_{T} (GeV/c)  (#eta = 0)"); fr->GetYaxis()->SetTitle("#sigma(p_{T})/p_{T}");
    fr->GetXaxis()->SetTitleSize(0.05); fr->GetYaxis()->SetTitleSize(0.05); fr->GetYaxis()->SetTitleOffset(1.5);
    TLegend* lg = new TLegend(0.19, 0.68, 0.96, 0.95); lg->SetBorderSize(0); lg->SetFillColor(kWhite); lg->SetFillStyle(1001); lg->SetTextSize(0.034);
    std::vector<TF1*> lines;
    for (auto& f : fits) {
      std::vector<double> bj = pr::best(Form("%s/logs/fitMomRes_%s.txt", note, f.tag), "Ks");
      TFile ff(Form("%s/fits/%s/fitMomRes_%s.root", note, f.tag, f.tag));
      TH2D* hj = (TH2D*)ff.Get("chi2_ab_Ks");
      if (bj.empty() || !hj) { printf("MISSING %s\n", f.tag); continue; }
      std::vector<double> lo(n, 1e9), hi(n, -1e9);
      for (int ia = 1; ia <= hj->GetNbinsX(); ++ia)
        for (int ib = 1; ib <= hj->GetNbinsY(); ++ib) {
          if (hj->GetBinContent(ia, ib) > 2.3) continue;
          const double a = hj->GetXaxis()->GetBinCenter(ia), b = hj->GetYaxis()->GetBinCenter(ib);
          for (int k = 0; k < n; ++k) {
            const double v = pr::res(a, b, 0.1 + 2.9 * k / (n - 1), m[i]);
            lo[k] = std::min(lo[k], v); hi[k] = std::max(hi[k], v);
          }
        }
      TGraph* band = new TGraph(2 * n);
      for (int k = 0; k < n; ++k) {
        const double pt = 0.1 + 2.9 * k / (n - 1);
        band->SetPoint(k, pt, hi[k]); band->SetPoint(2 * n - 1 - k, pt, lo[k]);
      }
      band->SetFillColorAlpha(f.col, 0.2); band->SetLineWidth(0); band->Draw("F");
      TF1* fb = new TF1(Form("fb_%s_%d", f.tag, i), "sqrt(([0]*x)^2 + ([1]*sqrt(x*x+[2]*[2])/x)^2)", 0.1, 3.0);
      fb->SetParameters(bj[0], bj[1], m[i]); fb->SetLineColor(f.col); fb->SetLineStyle(f.lsty); fb->SetLineWidth(3);
      lines.push_back(fb);
      lg->AddEntry(fb, Form("%s: a=%.2f%%/GeV, b=%.2f%%, s=%.3f", f.label, 100 * bj[0], 100 * bj[1], bj[2]), "l");
      printf("CURVE %s %s %.4f %.4f %.4f\n", f.tag, sp[i], pr::res(bj[0], bj[1], 1, m[i]), pr::res(bj[0], bj[1], 2, m[i]),
             pr::res(bj[0], bj[1], 3, m[i]));
      if (i == 0) {
        const int k = std::lround((2.0 - 0.1) / 2.9 * (n - 1));
        printf("BAND %s pi 2.0 %.4f %.4f\n", f.tag, lo[k], hi[k]);
      }
    }
    for (auto* l : lines) l->Draw("same");
    TH1* hs = (TH1*)fi.Get(Form("hMomRes_%s", sp[i]));
    if (hs) {
      TH1* h1 = (TH1*)hs->Clone(Form("h1_%d", i));
      h1->SetLineColor(kGray + 2); h1->SetLineWidth(2); h1->SetMarkerSize(0);
      h1->Draw("hist same ][");
      lg->AddEntry(h1, "full simulation (DDSim input)", "l");
    }
    lg->Draw();
    TLatex l; l.SetNDC(); l.SetTextSize(0.06); l.SetTextFont(42); l.DrawLatex(0.20, 0.60, lab[i]);
  }
  c.SaveAs(Form("%s/figures/fig_res.pdf", note));
  c.SaveAs(Form("%s/figures/fig_res.png", note));

  // ---- Delta chi2 maps of the Ks fits
  TCanvas c2("c2", "", 1300, 1050);
  c2.Divide(2, 2);
  for (int k = 0; k < 4; ++k) {
    c2.cd(k + 1);
    gPad->SetLeftMargin(0.15); gPad->SetRightMargin(0.15); gPad->SetBottomMargin(0.13); gPad->SetTopMargin(0.09);
    TFile ff(Form("%s/fits/%s/fitMomRes_%s.root", note, fits[k].tag, fits[k].tag));
    TH2D* h = (TH2D*)ff.Get("chi2_ab_Ks");
    if (!h) continue;
    h = (TH2D*)h->Clone(Form("m%d", k)); h->SetDirectory(nullptr);
    h->SetMaximum(25); h->SetMinimum(0);
    h->GetXaxis()->SetTitle("a (1/(GeV/c))"); h->GetYaxis()->SetTitle("b");
    h->GetXaxis()->SetTitleSize(0.05); h->GetYaxis()->SetTitleSize(0.05); h->GetYaxis()->SetTitleOffset(1.4);
    std::vector<double> b0 = pr::best(Form("%s/logs/fitMomRes_%s.txt", note, fits[k].tag), "Ks");
    h->GetXaxis()->SetRangeUser(b0[0] - 0.016, b0[0] + 0.016);
    h->GetYaxis()->SetRangeUser(std::max(0.004, b0[1] - 0.009), b0[1] + 0.009);
    h->GetZaxis()->SetTitle("#Delta#chi^{2} (min over s)");
    h->Draw("colz");
    std::vector<double> bk = pr::best(Form("%s/logs/fitMomRes_%s.txt", note, fits[k].tag), "Ks");
    TGraph* g = new TGraph(1); g->SetPoint(0, bk[0], bk[1]); g->SetMarkerStyle(29); g->SetMarkerSize(2.2); g->SetMarkerColor(kRed + 1);
    g->Draw("P same");
    TLatex l; l.SetNDC(); l.SetTextSize(0.05); l.SetTextFont(42);
    l.DrawLatex(0.16, 0.93, Form("K^{0}_{S}, %s: a=%.4f, b=%.3f, s=%.3f", fits[k].label, bk[0], bk[1], bk[2]));
  }
  c2.SaveAs(Form("%s/figures/fig_chi2.pdf", note));
  c2.SaveAs(Form("%s/figures/fig_chi2.png", note));
}
