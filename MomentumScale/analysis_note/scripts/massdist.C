// massdist.C -- invariant-mass distributions of Ks and Lambda+anti-Lambda before (recomputed, no kappa) and after the
// kappa pT correction, from the same QA-selected candidates. Also the KFParticle mass as reference.
//   root -l -b -q 'scripts/massdist.C("<note dir>","<merged V0mass dir>")'
// Figures (figures/, PNG only: ROOT PDF output renders the line widths inconsistently in this layout):
//   fig_massdist_int       2x2: Ks / Lambda (rows) x period A / B (columns), pT-integrated (Ks 0.2-4.0, Lambda 0.5-4.0 GeV/c)
//   fig_massdist_Ks_<P>    Ks per V0-pT bin of the momentum fit, period P
//   fig_massdist_La_<P>    Lambda+anti-Lambda per V0-pT bin, period P
// Binning: Ks 2 MeV, Lambda 1 MeV (pT-integrated). Same candidates in all versions; the V0 pT of the binning is that version's own pT.
// Prints MDIST rows: channel period version ptlo pthi entries(in the plotted mass range) bin width (MeV).
// Per-pT-bin panels with few entries are rebinned (Ks 4 or 8 MeV, Lambda 2 or 4 MeV; same width for all versions of a panel).
#include <TFile.h>
#include <TH2F.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TLine.h>
#include <TStyle.h>
#include <cstdio>
#include <vector>

namespace md {
TH1D* proj(TFile& f, bool ks, double plo, double phi, const char* nm, int extra)
{
  TH2F* h = (TH2F*)f.Get(ks ? "hKs_qa" : "hLa_qa")->Clone(Form("%s_2d", nm));
  if (!ks) h->Add((TH2F*)f.Get("hLaBar_qa"));
  TH1D* p = h->ProjectionY(nm, h->GetXaxis()->FindBin(plo + 1e-6), h->GetXaxis()->FindBin(phi - 1e-6));
  p->SetDirectory(nullptr);
  p->Rebin((ks ? 2 : 4) * extra);
  p->GetXaxis()->SetRangeUser(ks ? 0.44 : 1.09, ks ? 0.56 : 1.145);
  return p;
}
void style(TH1D* h, int col, int ls, int lw) { h->SetLineColor(col); h->SetLineStyle(ls); h->SetLineWidth(lw); h->SetMarkerSize(0); }
}  // namespace md

void massdist(const char* note, const char* dataDir)
{
  gStyle->SetOptStat(0); gStyle->SetOptTitle(0);
  gStyle->SetTitleSize(0.055, "XY"); gStyle->SetLabelSize(0.048, "XY"); gStyle->SetTitleOffset(1.35, "Y");
  const char* per[2] = {"A", "B"};
  const char* ver[3] = {"kfp", "raw", "kappa"};
  const char* vlab[3] = {"KFParticle mass", "recomputed, no #kappa", "recomputed, with #kappa"};
  const int col[3] = {kGray + 2, kBlue + 1, kRed + 1}, lst[3] = {2, 1, 1}, lwd[3] = {2, 3, 3};
  const std::vector<double> ptK = {0.2, 0.5, 0.7, 0.9, 1.2, 1.6, 2.2, 4.0}, ptL = {0.5, 0.9, 1.2, 1.6, 2.4, 4.0};
  TFile* f[2][3];
  for (int p = 0; p < 2; ++p) for (int v = 0; v < 3; ++v) f[p][v] = TFile::Open(Form("%s/V0mass_%s_%s.root", dataDir, per[p], ver[v]));

  TH1D* hleg[3] = {nullptr, nullptr, nullptr}; TLine* lleg = nullptr;
  auto panel = [&](bool ks, int p, double lo, double hi, const char* tag, bool legend, double tsz, bool adapt) {
    TH1D* h[3]; double ymax = 0;
    int extra = 1;
    if (adapt) {
      TH1D* t = md::proj(*f[p][1], ks, lo, hi, Form("t_%s_%s_%s", tag, per[p], ks ? "K" : "L"), 1);
      const double n = t->Integral();
      extra = n < 1500 ? 4 : (n < 4000 ? 2 : 1);
      delete t;
    }
    const double bw = (ks ? 2.0 : 1.0) * extra;
    for (int v = 0; v < 3; ++v) {
      h[v] = md::proj(*f[p][v], ks, lo, hi, Form("m_%s_%s_%d_%s", tag, per[p], v, ks ? "K" : "L"), extra);
      md::style(h[v], col[v], lst[v], lwd[v]);
      ymax = std::max(ymax, h[v]->GetMaximum());
      printf("MDIST %s %s %s %.1f %.1f %.0f %.0f\n", ks ? "Ks" : "Lambda", per[p], ver[v], lo, hi, h[v]->Integral(), bw);
      hleg[v] = h[v];
    }
    gPad->SetLeftMargin(0.17); gPad->SetRightMargin(0.03); gPad->SetTopMargin(legend ? 0.08 : 0.17); gPad->SetBottomMargin(0.15);
    h[1]->SetMaximum((legend ? 1.30 : 1.10) * ymax); h[1]->SetMinimum(0);
    h[1]->SetXTitle(ks ? "M(#pi^{+}#pi^{-}) (GeV/c^{2})" : "M(p#pi) (GeV/c^{2})");
    h[1]->SetYTitle(Form("candidates / %.0f MeV/c^{2}", bw));
    h[1]->GetXaxis()->SetNdivisions(505);
    h[1]->GetXaxis()->SetTitleSize(tsz); h[1]->GetYaxis()->SetTitleSize(tsz);
    h[1]->GetXaxis()->SetLabelSize(0.9 * tsz); h[1]->GetYaxis()->SetLabelSize(0.9 * tsz);
    h[1]->Draw("hist"); h[0]->Draw("hist same"); h[1]->Draw("hist same"); h[2]->Draw("hist same");
    const double m0 = ks ? 0.497611 : 1.115683;
    TLine* l = new TLine(m0, 0, m0, (legend ? 1.30 : 1.10) * ymax); l->SetLineStyle(3); l->SetLineColor(kBlack); l->SetLineWidth(2); l->Draw(); lleg = l;
    TLatex t; t.SetNDC(); t.SetTextSize(tsz);
    if (legend) {
      t.DrawLatex(0.20, 0.85, Form("%s, period %s", ks ? "K^{0}_{S}" : "#Lambda+#bar{#Lambda}", per[p]));
      t.DrawLatex(0.20, 0.78, Form("%.1f < p_{T} < %.1f GeV/c", lo, hi));
    } else {   // per-bin panels: label above the frame
      t.SetTextSize(0.9 * tsz); t.DrawLatex(0.17, 0.92, Form("%s (%s):  %.1f < p_{T} < %.1f GeV/c", ks ? "K^{0}_{S}" : "#Lambda+#bar{#Lambda}", per[p], lo, hi));
    }
    if (legend) {
      TLegend* lg = new TLegend(0.60, 0.70, 0.97, 0.91); lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.85 * tsz);
      for (int v : {1, 2, 0}) lg->AddEntry(h[v], vlab[v], "l");
      lg->AddEntry(l, "PDG", "l"); lg->Draw();
    }
  };

  TCanvas c1("c1", "", 2000, 1500); c1.Divide(2, 2);
  for (int c = 0; c < 2; ++c) for (int p = 0; p < 2; ++p) {
    c1.cd(2 * c + p + 1);
    panel(c == 0, p, c == 0 ? 0.2 : 0.5, 4.0, "int", true, 0.050, false);
  }
  c1.SaveAs(Form("%s/figures/fig_massdist_int.png", note));

  for (int c = 0; c < 2; ++c) for (int p = 0; p < 2; ++p) {
    const bool ks = c == 0; const auto& pt = ks ? ptK : ptL;
    const int nb = pt.size() - 1, nx = ks ? 4 : 3, ny = 2;
    TCanvas cv("cv", "", 2400, 1250); cv.Divide(nx, ny);
    for (int j = 0; j < nb; ++j) { cv.cd(j + 1); panel(ks, p, pt[j], pt[j + 1], Form("b%d", j), false, 0.062, true); }
    cv.cd(nb + 1);
    TLegend* lg = new TLegend(0.05, 0.25, 0.95, 0.85); lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.085);
    lg->SetHeader(Form("%s, period %s", ks ? "K^{0}_{S}" : "#Lambda+#bar{#Lambda}", per[p]));
    for (int v : {1, 2, 0}) lg->AddEntry(hleg[v], vlab[v], "l");
    lg->AddEntry(lleg, "PDG mass", "l"); lg->Draw();
    cv.SaveAs(Form("%s/figures/fig_massdist_%s_%s.png", note, ks ? "Ks" : "La", per[p]));
  }
}
