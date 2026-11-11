// plotKappaCompare.C -- before/after the kappa pT re-calibration (run 79509 seg 0-30, loosecov):
// (1) Gaussian peak position and width vs V0 pT (data), (2) fitted sigma(pT)/pT for pions and protons from the
// (a, b, s) template fits, (3) pT-integrated mass peaks. Fit results are passed as arguments (from FIT lines).
#include <TFile.h>
#include <TH1D.h>
#include <TH2F.h>
#include <TF1.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TLine.h>
#include <TStyle.h>
#include <TPad.h>
#include <cmath>
#include <vector>

namespace pk {
TGraphErrors* peak(const char* file, const char* h1, const char* h2, const std::vector<double>& pt, double m0, double w,
                   double lo, double hi, bool width, int rebin)
{
  TFile F(file);
  TH2F* h = (TH2F*)F.Get(h1)->Clone(); h->SetDirectory(nullptr);
  if (h2) h->Add((TH2F*)F.Get(h2));
  TGraphErrors* g = new TGraphErrors();
  for (size_t j = 0; j + 1 < pt.size(); ++j) {
    TH1D* p = h->ProjectionY(Form("pp_%s_%zu_%d", h1, j, width), h->GetXaxis()->FindBin(pt[j] + 1e-6), h->GetXaxis()->FindBin(pt[j + 1] - 1e-6));
    p->Rebin(rebin);
    if (p->Integral() < 100) continue;
    TF1 f("f", "gaus(0)+pol1(3)", lo, hi);
    f.SetParameters(p->GetMaximum(), m0, w, 0, 0); f.SetParLimits(2, 0.0003, 0.02); f.SetParLimits(1, m0 - 0.01, m0 + 0.01);
    p->Fit(&f, "QNR");
    const int n = g->GetN();
    const double v = width ? std::fabs(f.GetParameter(2)) : f.GetParameter(1), e = width ? f.GetParError(2) : f.GetParError(1);
    if (e > (width ? 0.004 : 0.004)) continue;   // drop unconstrained fits
    g->SetPoint(n, 0.5 * (pt[j] + pt[j + 1]), 1e3 * v);
    g->SetPointError(n, 0.5 * (pt[j + 1] - pt[j]), 1e3 * e);
  }
  return g;
}
}  // namespace pk

void plotKappaCompare(double aB, double bB, double sB, double aA, double bA, double sA, const char* out = "output/momfit_kappa/kappa_compare")
{
  gStyle->SetOptStat(0);
  gStyle->SetTitleSize(0.05, "XY"); gStyle->SetLabelSize(0.045, "XY"); gStyle->SetTitleOffset(1.7, "Y"); gStyle->SetTitleOffset(1.2, "X");
  const char* fB = "data/kappa/V0mass_79509seg000-030_mode1.root";
  const char* fA = "data/kappa/V0mass_79509seg000-030_mode2.root";
  std::vector<double> ptK = {0.2, 0.5, 0.7, 0.9, 1.2, 1.6, 2.2, 4.0}, ptL = {0.5, 0.9, 1.2, 1.6, 2.4, 4.0};
  TCanvas c("c", "", 1500, 1000);
  c.Divide(3, 2);
  struct P { const char* t; const char* h1; const char* h2; std::vector<double>* pt; double m0, w, lo, hi; bool width; double y0, y1; int rebin; };
  P ps[4] = {{"K^{0}_{S} peak position", "hKs_qa", nullptr, &ptK, 0.4976, 0.005, 0.46, 0.54, false, 490, 502, 2},
             {"K^{0}_{S} Gaussian width", "hKs_qa", nullptr, &ptK, 0.4976, 0.005, 0.46, 0.54, true, 0, 12, 2},
             {"#Lambda+#bar{#Lambda} peak position", "hLa_qa", "hLaBar_qa", &ptL, 1.1157, 0.002, 1.103, 1.130, false, 1112, 1119, 4},
             {"#Lambda+#bar{#Lambda} Gaussian width", "hLa_qa", "hLaBar_qa", &ptL, 1.1157, 0.002, 1.103, 1.130, true, 0, 5, 4}};
  int pad = 1;
  for (auto& p : ps) {
    c.cd(pad == 3 ? 4 : pad == 4 ? 5 : pad); ++pad;
    gPad->SetLeftMargin(0.19); gPad->SetTopMargin(0.12); gPad->SetBottomMargin(0.15); gPad->SetRightMargin(0.04);
    TH1F* fr = gPad->DrawFrame(0, p.y0, 4.0, p.y1);
    fr->SetXTitle("V0 p_{T} (GeV/c)"); fr->SetYTitle(p.width ? "#sigma (MeV/c^{2})" : "M (MeV/c^{2})");
    TGraphErrors* gB = pk::peak(fB, p.h1, p.h2, *p.pt, p.m0, p.w, p.lo, p.hi, p.width, p.rebin);
    TGraphErrors* gA = pk::peak(fA, p.h1, p.h2, *p.pt, p.m0, p.w, p.lo, p.hi, p.width, p.rebin);
    gB->SetMarkerStyle(24); gB->SetMarkerColor(kBlue + 1); gB->SetLineColor(kBlue + 1);
    gA->SetMarkerStyle(20); gA->SetMarkerColor(kRed + 1); gA->SetLineColor(kRed + 1);
    if (!p.width) { TLine* l = new TLine(0, 1e3 * p.m0, 4.0, 1e3 * p.m0); l->SetLineStyle(2); l->Draw(); }
    gB->Draw("P same"); gA->Draw("P same");
    TLatex t; t.SetNDC(); t.SetTextSize(0.05); t.DrawLatex(0.15, 0.91, p.t);
    TLegend* lg = new TLegend(0.45, 0.72, 0.88, 0.86); lg->SetBorderSize(0); lg->SetTextSize(0.04);
    lg->AddEntry(gB, "before (KFP momenta)", "p"); lg->AddEntry(gA, "after #kappa correction", "p"); lg->Draw();
  }
  // fitted single-track resolution
  c.cd(3); gPad->SetLeftMargin(0.19); gPad->SetTopMargin(0.12); gPad->SetBottomMargin(0.15); gPad->SetRightMargin(0.04);
  TH1F* fr = gPad->DrawFrame(0, 0, 4, 18);
  fr->SetXTitle("track p_{T} (GeV/c)"); fr->SetYTitle("#sigma_{p_{T}}/p_{T} (%)");
  auto res = [](double a, double b, double m, const char* n) {
    TF1* f = new TF1(n, Form("100*sqrt(pow(%g*x,2)+pow(%g,2)*(1+%g*%g/(x*x)))", a, b, m, m), 0.1, 4);   // p ~ pT at eta 0
    return f;
  };
  TF1* r[4] = {res(aB, bB, 0.1396, "rBpi"), res(aA, bA, 0.1396, "rApi"), res(aB, bB, 0.9383, "rBp"), res(aA, bA, 0.9383, "rAp")};
  int col[4] = {kBlue + 1, kRed + 1, kBlue + 1, kRed + 1}, sty[4] = {1, 1, 2, 2};
  TLegend* lg = new TLegend(0.21, 0.66, 0.95, 0.87); lg->SetBorderSize(0); lg->SetTextSize(0.033);
  const char* lab[4] = {Form("#pi before: a=%.2f%%/GeV b=%.2f%% s=%.3f", 100 * aB, 100 * bB, sB), Form("#pi after:   a=%.2f%%/GeV b=%.2f%% s=%.3f", 100 * aA, 100 * bA, sA), "p before", "p after"};
  for (int i = 0; i < 4; ++i) { r[i]->SetLineColor(col[i]); r[i]->SetLineStyle(sty[i]); r[i]->Draw("same"); lg->AddEntry(r[i], lab[i], "l"); }
  lg->Draw();
  TLatex t; t.SetNDC(); t.SetTextSize(0.045);
  t.DrawLatex(0.18, 0.91, "fitted #sigma_{p_{T}}/p_{T} (#eta = 0)");
  // integrated peaks
  for (int k = 0; k < 2; ++k) {
    c.cd(k == 0 ? 6 : 6);
  }
  c.cd(6); gPad->SetLeftMargin(0.19); gPad->SetTopMargin(0.12); gPad->SetBottomMargin(0.15); gPad->SetRightMargin(0.04);
  {
    TFile FB(fB), FA(fA);
    TH1D* hb = ((TH2F*)FB.Get("hKs_qa"))->ProjectionY("ksB", 5, 80); hb->SetDirectory(nullptr);
    TH1D* ha = ((TH2F*)FA.Get("hKs_qa"))->ProjectionY("ksA", 5, 80); ha->SetDirectory(nullptr);
    hb->Rebin(2); ha->Rebin(2);
    hb->GetXaxis()->SetRangeUser(0.45, 0.55); hb->SetLineColor(kBlue + 1); ha->SetLineColor(kRed + 1);
    hb->SetMaximum(1.15 * std::max(hb->GetMaximum(), ha->GetMaximum())); hb->SetMinimum(0);
    hb->SetXTitle("M_{#pi#pi} (GeV/c^{2})"); hb->GetXaxis()->SetNdivisions(505); hb->UseCurrentStyle(); hb->SetLineColor(kBlue + 1); hb->SetTitle("");
    hb->Draw("hist"); ha->Draw("hist same");
    TLine* l = new TLine(0.49761, 0, 0.49761, hb->GetMaximum()); l->SetLineStyle(2); l->Draw();
    TLatex tt; tt.SetNDC(); tt.SetTextSize(0.05); tt.DrawLatex(0.15, 0.91, "K^{0}_{S}, 0.2<p_{T}<4 GeV/c, QA cuts");
    TLegend* lg2 = new TLegend(0.55, 0.72, 0.88, 0.86); lg2->SetBorderSize(0); lg2->SetTextSize(0.04);
    lg2->AddEntry(hb, "before", "l"); lg2->AddEntry(ha, "after #kappa", "l"); lg2->Draw();
  }
  c.SaveAs(Form("%s.png", out));
  c.SaveAs(Form("%s.pdf", out));
}
