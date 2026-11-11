// d0_bkgcheck.C : for one selection, fit a fixed-width gaussian at many trial masses with pol2 / pol3 / pol4 background
// A real D0 must stand out at 1.865 and be stable against the background model.
#include <ROOT/RDataFrame.hxx>
#include <TCanvas.h>
#include <TF1.h>
#include <TGraph.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TLine.h>
#include <TStyle.h>

#include <cmath>
#include <cstdio>

void d0_bkgcheck(const char *cut = "D0_DIRA>0.995 && D0_decayLength/D0_decayLengthErr>2", double sig = 0.015)
{
  ROOT::EnableImplicitMT(6);
  gStyle->SetOptStat(0);
  ROOT::RDataFrame d("DecayTree", "../ntuples/d0skim/D0skim_prod10pct_*.root");
  auto hp = d.Filter(cut).Histo1D({"h", Form("%s;m(K#pi) (GeV/c^{2});counts / 5 MeV", cut), 100, 1.6, 2.1}, "D0_mass");
  TH1D *h = hp.GetPtr();
  TCanvas c("c", "c", 1400, 600);
  c.Divide(2, 1);
  TGraph g[3];
  int col[3] = {kBlack, kBlue + 1, kRed + 1};
  for (int ord = 2; ord <= 4; ++ord)
  {
    int ip = 0;
    for (double m = 1.74; m <= 1.99; m += 0.005)
    {
      TF1 f("f", Form("[0]*exp(-0.5*((x-[1])/[2])^2)+pol%d(3)", ord), 1.66, 2.08);
      f.SetParameter(0, 0);
      f.FixParameter(1, m);
      f.FixParameter(2, sig);
      for (int k = 3; k < 4 + ord; ++k) f.SetParameter(k, k == 3 ? h->GetBinContent(50) : 0);
      h->Fit(&f, "QN0RL");
      h->Fit(&f, "QN0RL");
      double z = f.GetParError(0) > 0 ? f.GetParameter(0) / f.GetParError(0) : 0;
      g[ord - 2].SetPoint(ip++, m, z);
      if (std::fabs(m - 1.865) < 0.0026)
        printf("pol%d: at m=%.3f S/err = %.2f  chi2/ndf = %.2f\n", ord, m, z, f.GetChisquare() / f.GetNDF());
    }
  }
  // print the pol3 local-significance profile
  for (int i = 0; i < g[1].GetN(); i += 2) printf("m=%.3f  pol2 %.2f  pol3 %.2f  pol4 %.2f\n", g[0].GetX()[i], g[0].GetY()[i], g[1].GetY()[i], g[2].GetY()[i]);
  c.cd(1);
  h->SetMinimum(0);
  h->Draw("e");
  c.cd(2);
  TLegend l(0.15, 0.75, 0.45, 0.9);
  for (int k = 0; k < 3; ++k)
  {
    g[k].SetLineColor(col[k]);
    g[k].SetLineWidth(2);
    g[k].Draw(k ? "L" : "AL");
    l.AddEntry(&g[k], Form("pol%d background", k + 2), "l");
  }
  g[0].SetTitle(Form("local S/#sigma_{S} vs trial mass (gaus #sigma=%.0f MeV);trial mass (GeV/c^{2});S/#sigma_{S}", sig * 1000));
  g[0].GetYaxis()->SetRangeUser(-5, 6);
  TLine ln(1.8648, -5, 1.8648, 6);
  ln.SetLineStyle(2);
  ln.Draw();
  l.Draw();
  c.SaveAs("d0_bkgcheck.png");
}
