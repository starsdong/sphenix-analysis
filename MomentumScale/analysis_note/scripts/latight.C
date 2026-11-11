// latight.C -- Lambda+anti-Lambda mass with the QA and with the tight selection, recomputed with kappa, pT 0.5-4 GeV/c,
// periods A and B (condor 481293, data/kappa_v2). Draws figures/fig_latight_mass.png.
//   root -l -b -q 'scripts/latight.C("<note dir>","<dir with V0mass_{A,B}_kappa.root>")'
#include <TFile.h>
#include <TH2.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TLine.h>
#include <TStyle.h>

void latight(const char* note, const char* dataDir)
{
  gStyle->SetOptStat(0); gStyle->SetOptTitle(0);
  TCanvas cv("cv", "", 2000, 800); cv.Divide(2, 1);
  const char* per[2] = {"A", "B"};
  for (int p = 0; p < 2; ++p) {
    TFile f(Form("%s/V0mass_%s_kappa.root", dataDir, per[p]));
    TH1D* h[2];
    for (int c = 0; c < 2; ++c) {
      TH2* x = (TH2*)f.Get(c ? "hLa_tight" : "hLa_qa")->Clone(Form("x%d%d", p, c));
      x->Add((TH2*)f.Get(c ? "hLaBar_tight" : "hLaBar_qa"));
      h[c] = x->ProjectionY(Form("h%d%d", p, c), x->GetXaxis()->FindBin(0.5001), x->GetXaxis()->FindBin(3.9999));
      h[c]->SetDirectory(nullptr); h[c]->Rebin(4);
      h[c]->GetXaxis()->SetRangeUser(1.09, 1.145);
      h[c]->SetLineColor(c ? kRed + 1 : kBlue + 1); h[c]->SetLineWidth(3);
    }
    cv.cd(p + 1);
    gPad->SetLeftMargin(0.17); gPad->SetRightMargin(0.03); gPad->SetTopMargin(0.07); gPad->SetBottomMargin(0.13);
    h[0]->SetMaximum(1.35 * h[0]->GetMaximum()); h[0]->SetMinimum(0);
    h[0]->SetXTitle("M(p#pi) (GeV/c^{2})"); h[0]->SetYTitle("candidates / 1 MeV/c^{2}");
    for (TH1D* x : h) { x->GetXaxis()->SetTitleSize(0.05); x->GetYaxis()->SetTitleSize(0.05); x->GetXaxis()->SetLabelSize(0.045); x->GetYaxis()->SetLabelSize(0.045); x->GetYaxis()->SetTitleOffset(1.5); x->GetXaxis()->SetNdivisions(505); }
    h[0]->Draw("hist"); h[1]->Draw("hist same");
    TLine* l = new TLine(1.115683, 0, 1.115683, h[0]->GetMaximum()); l->SetLineStyle(3); l->SetLineWidth(2); l->Draw();
    TLatex t; t.SetNDC(); t.SetTextSize(0.05);
    t.DrawLatex(0.21, 0.86, Form("#Lambda+#bar{#Lambda}, period %s, with #kappa", per[p]));
    t.DrawLatex(0.21, 0.79, "0.5 < p_{T} < 4.0 GeV/c");
    TLegend* lg = new TLegend(0.62, 0.72, 0.97, 0.91); lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.042);
    lg->AddEntry(h[0], "QA selection", "l"); lg->AddEntry(h[1], "tight selection", "l"); lg->AddEntry(l, "PDG mass", "l"); lg->Draw();
  }
  cv.SaveAs(Form("%s/figures/fig_latight_mass.png", note));
}
