// dedx_map.C : dE/dx vs p for D0-skim daughters, and basic distributions
#include <ROOT/RDataFrame.hxx>
#include <TCanvas.h>
#include <TH2D.h>
#include <TStyle.h>
#include <TLatex.h>
#include <TSystem.h>

void dedx_map()
{
  ROOT::EnableImplicitMT(6);
  gStyle->SetOptStat(0);
  ROOT::RDataFrame d("DecayTree", "../ntuples/d0skim/D0skim_prod10pct_*.root");
  auto d2 = d.Define("p1", "track_1_pT*cosh(track_1_pseudorapidity)")
                .Define("p2", "track_2_pT*cosh(track_2_pseudorapidity)");
  auto h1 = d2.Histo2D({"h1", "track 1 (#pi hyp.);p (GeV/c);dE/dx (a.u.)", 150, 0, 3, 200, 0, 1500}, "p1", "track_1_calculated_dEdx");
  auto h2 = d2.Histo2D({"h2", "track 2 (K hyp.);p (GeV/c);dE/dx (a.u.)", 150, 0, 3, 200, 0, 1500}, "p2", "track_2_calculated_dEdx");
  auto hm = d2.Histo1D({"hm", "all skim;m(K#pi) (GeV/c^{2});", 100, 1.5, 2.2}, "D0_mass");
  auto hpt = d2.Histo1D({"hpt", ";D0 p_{T};", 100, 0, 8}, "D0_pT");
  auto hls = d2.Define("ls", "D0_decayLength/D0_decayLengthErr").Histo1D({"hls", ";L/#sigma_{L};", 100, -5, 20}, "ls");
  auto hdira = d2.Histo1D({"hdira", ";DIRA;", 100, 0.95, 1.0}, "D0_DIRA");
  TCanvas c("c", "c", 1500, 900);
  c.Divide(3, 2);
  c.cd(1); gPad->SetLogz(); h1->DrawClone("colz");
  c.cd(2); gPad->SetLogz(); h2->DrawClone("colz");
  c.cd(3); hm->DrawClone();
  c.cd(4); gPad->SetLogy(); hpt->DrawClone();
  c.cd(5); gPad->SetLogy(); hls->DrawClone();
  c.cd(6); gPad->SetLogy(); hdira->DrawClone();
  c.SaveAs("dedx_map.png");
  // print dE/dx quantiles in p slices (to locate pion / kaon / proton bands)
  for (double plo : {0.5, 0.6, 0.7, 0.8, 1.0, 1.5})
  {
    auto *py = h1->ProjectionY("py", h1->GetXaxis()->FindBin(plo + 1e-3), h1->GetXaxis()->FindBin(plo + 0.1 - 1e-3));
    double q[5], pr[5] = {0.05, 0.25, 0.5, 0.75, 0.95};
    py->GetQuantiles(5, q, pr);
    printf("p %.1f-%.1f: N=%.0f dEdx quantiles 5/25/50/75/95%% = %.0f %.0f %.0f %.0f %.0f\n", plo, plo + 0.1, py->GetEntries(), q[0], q[1], q[2], q[3], q[4]);
  }
  printf("mass range %.3f-%.3f, entries %.0f\n", hm->GetXaxis()->GetXmin(), hm->GetXaxis()->GetXmax(), hm->GetEntries());
}
