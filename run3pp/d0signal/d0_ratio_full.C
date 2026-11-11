// d0_ratio.C : mass spectrum of pointing candidates divided by non-pointing (DIRA 0.95-0.98, background-dominated)
// control sample with the same daughter kinematics. A D0 signal shows up as a peak in the ratio.
#include <ROOT/RDataFrame.hxx>
#include <TCanvas.h>
#include <TF1.h>
#include <TH1D.h>
#include <TLatex.h>
#include <TLine.h>
#include <TStyle.h>
#include <TSystem.h>

#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

void d0_ratio_full()
{
  ROOT::EnableImplicitMT(6);
  gStyle->SetOptStat(0);
  ROOT::RDataFrame d0("DecayTree", "../ntuples/d0tight_full/D0tight_pT1p5_run*.root");
  auto d = d0.Define("ls", "D0_decayLength/D0_decayLengthErr");
  std::string ctrl = "D0_DIRA>0.95 && D0_DIRA<0.98";
  struct S
  {
    std::string name, sig, ctl;
  };
  std::vector<S> sets = {
      {"pT>1.5: DIRA>0.995, L/#sigma>2", "D0_DIRA>0.995 && ls>2", ctrl},
      {"pT>1.5: DIRA>0.998, L/#sigma>2", "D0_DIRA>0.998 && ls>2", ctrl},
      {"pT>1.5: DIRA>0.995, L/#sigma>4", "D0_DIRA>0.995 && ls>4", ctrl},
      {"pT 1.5-2.5: DIRA>0.995, L/#sigma>2", "D0_DIRA>0.995 && ls>2 && D0_pT<2.5", ctrl + " && D0_pT<2.5"},
      {"pT>2.5: DIRA>0.995, L/#sigma>2", "D0_DIRA>0.995 && ls>2 && D0_pT>2.5", ctrl + " && D0_pT>2.5"},
      {"control: DIRA 0.98-0.99", "D0_DIRA>0.98 && D0_DIRA<0.99", ctrl},
  };
  std::vector<ROOT::RDF::RResultPtr<TH1D>> hs, hc;
  for (size_t i = 0; i < sets.size(); ++i)
  {
    hs.push_back(d.Filter(sets[i].sig).Histo1D({Form("s%zu", i), "", 50, 1.6, 2.1}, "D0_mass"));
    hc.push_back(d.Filter(sets[i].ctl).Histo1D({Form("c%zu", i), "", 50, 1.6, 2.1}, "D0_mass"));
  }
  TCanvas c("c", "c", 1500, 900);
  c.Divide(3, 2);
  for (size_t i = 0; i < sets.size(); ++i)
  {
    auto *r = new TH1D(Form("r%zu", i), "", 50, 1.6, 2.1);
    for (int b = 1; b <= 50; ++b)
    {
      double ns = hs[i]->GetBinContent(b), nc = hc[i]->GetBinContent(b);
      if (ns <= 0 || nc <= 0) continue;
      r->SetBinContent(b, ns / nc);
      r->SetBinError(b, ns / nc * std::sqrt(1. / ns + 1. / nc));
    }
    r->Scale(1. / r->Integral() * r->GetNbinsX());
    r->SetTitle((sets[i].name + ";m(K#pi) (GeV/c^{2});ratio (normalised)").c_str());
    // fit: gaus (free) + pol2
    TF1 f("f", "gaus(0)+pol2(3)", 1.65, 2.08);
    f.SetParameters(0.05, 1.865, 0.03, 1, 0, 0);
    f.SetParLimits(1, 1.82, 1.92);
    f.SetParLimits(2, 0.008, 0.08);
    r->Fit(&f, "QR0");
    TF1 f0("f0", "pol2", 1.65, 2.08);
    r->Fit(&f0, "QR0");
    double dchi2 = f0.GetChisquare() - f.GetChisquare();
    printf("%-58s Nsig %8.0f Nctl %8.0f | gaus amp %.4f+-%.4f mean %.4f+-%.4f sigma %.4f+-%.4f | chi2/ndf %.1f/%d (pol2 only %.1f/%d) dchi2 %.1f\n",
           sets[i].name.c_str(), hs[i]->GetEntries(), hc[i]->GetEntries(), f.GetParameter(0), f.GetParError(0), f.GetParameter(1), f.GetParError(1),
           f.GetParameter(2), f.GetParError(2), f.GetChisquare(), f.GetNDF(), f0.GetChisquare(), f0.GetNDF(), dchi2);
    c.cd(i + 1);
    r->Draw("e");
    f.SetLineColor(kRed);
    f.DrawCopy("same");
    f0.SetLineColor(kGray + 2);
    f0.SetLineStyle(2);
    f0.DrawCopy("same");
    TLine l(1.8648, r->GetMinimum(), 1.8648, r->GetMaximum());
    l.SetLineStyle(3);
    l.DrawClone();
    TLatex t;
    t.SetNDC();
    t.SetTextSize(0.045);
    t.DrawLatex(0.15, 0.85, Form("#mu=%.3f #sigma=%.0f MeV, #Delta#chi^{2}=%.0f", f.GetParameter(1), 1000 * f.GetParameter(2), dchi2));
  }
  c.SaveAs(gSystem->Getenv("RATIO_PNG") ? gSystem->Getenv("RATIO_PNG") : "d0_ratio.png");
}
