// d0_grid.C : grid scan of D0 selections; peak fixed at the PDG mass, sigma fixed (default 15 MeV); pol2 background
// usage: root -b -q 'd0_grid.C+(0.015)'     output: d0_grid_<sigma>.txt (sorted by S/err), d0_grid_best.png
#include <Math/Boost.h>
#include <Math/Vector4D.h>
#include <ROOT/RDataFrame.hxx>
#include <TCanvas.h>
#include <TF1.h>
#include <TH1D.h>
#include <TLatex.h>
#include <TStyle.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

void d0_grid(double sig = 0.015)
{
  ROOT::EnableImplicitMT(6);
  gStyle->SetOptStat(0);
  ROOT::RDataFrame d0("DecayTree", "../ntuples/d0skim/D0skim_prod10pct_*.root");
  // pion band (from d0_scan.C): median / width vs p, linear interpolation of the measured points
  auto med = [](float p) { p = std::min(std::max(p, 0.55f), 2.2f); return p < 0.9 ? 198 + (p - 0.55) / 0.35 * 40 : (p < 1.6 ? 238 + (p - 0.9) / 0.7 * 43 : 281 + (p - 1.6) / 0.6 * 11); };
  auto wid = [](float p) { p = std::min(std::max(p, 0.55f), 2.2f); return p < 0.9 ? 43 + (p - 0.55) / 0.35 * 9 : (p < 1.6 ? 52 + (p - 0.9) / 0.7 * 11 : 63 + (p - 1.6) / 0.6 * 5); };
  auto d = d0.Define("pK", "track_2_pT*cosh(track_2_pseudorapidity)")
               .Define("pP", "track_1_pT*cosh(track_1_pseudorapidity)")
               .Define("zK", [=](float p, float de) { return (de - med(p)) / wid(p); }, {"pK", "track_2_calculated_dEdx"})
               .Define("zP", [=](float p, float de) { return (de - med(p)) / wid(p); }, {"pP", "track_1_calculated_dEdx"})
               .Define("ls", "D0_decayLength/D0_decayLengthErr")
               .Define("cts",
                       [](float pt1, float eta1, float phi1, float pt2, float eta2, float phi2) {
                         using V = ROOT::Math::PtEtaPhiMVector;
                         V pi(pt1, eta1, phi1, 0.139570), k(pt2, eta2, phi2, 0.493677);
                         auto D = k + pi;
                         auto kr = ROOT::Math::Boost(D.BoostToCM())(k);
                         return (double) kr.Vect().Unit().Dot(D.Vect().Unit());
                       },
                       {"track_1_pT", "track_1_pseudorapidity", "track_1_phi", "track_2_pT", "track_2_pseudorapidity", "track_2_phi"});
  struct Cfg
  {
    std::string name, cut;
  };
  std::vector<Cfg> cfgs;
  for (double dira : {0.99, 0.995, 0.998, 0.999})
    for (double ls : {2., 3., 4., 5., 7.})
      for (double pt : {0., 1., 2., 3.})
        for (int pid : {0, 1})
          for (int cts : {0, 1})
          {
            std::string c = Form("D0_DIRA>%g && ls>%g && D0_pT>%g", dira, ls, pt);
            if (pid) c += " && (pK>0.9 || zK>2) && abs(zP)<3";
            if (cts) c += " && abs(cts)<0.8";
            cfgs.push_back({Form("DIRA>%.3f L/sL>%.0f pT>%.0f %s %s", dira, ls, pt, pid ? "PID" : "   ", cts ? "cts" : "   "), c});
          }
  std::vector<ROOT::RDF::RResultPtr<TH1D>> hs;
  for (size_t i = 0; i < cfgs.size(); ++i)
    hs.push_back(d.Filter(cfgs[i].cut).Histo1D({Form("g%zu", i), (cfgs[i].name + ";m(K#pi) (GeV/c^{2});counts / 5 MeV").c_str(), 100, 1.6, 2.1}, "D0_mass"));

  struct R
  {
    size_t i;
    double S, eS, B, N;
  };
  std::vector<R> res;
  const double mD0 = 1.86484;
  for (size_t i = 0; i < cfgs.size(); ++i)
  {
    TH1D *h = hs[i].GetPtr();
    if (h->Integral() < 100) continue;
    TF1 f(Form("f%zu", i), "[0]*exp(-0.5*((x-[1])/[2])^2)+pol2(3)", 1.70, 2.05);
    f.SetParameters(0, mD0, sig, h->GetBinContent(30), 0, 0);
    f.FixParameter(1, mD0);
    f.FixParameter(2, sig);
    h->Fit(&f, "QN0RL");  // amplitude allowed to go negative: unbiased significance
    double bw = h->GetBinWidth(1);
    double S = f.GetParameter(0) * sig * std::sqrt(2 * M_PI) / bw, eS = f.GetParError(0) * sig * std::sqrt(2 * M_PI) / bw;
    TF1 b("b", "pol2", 1.7, 2.05);
    b.SetParameters(f.GetParameter(3), f.GetParameter(4), f.GetParameter(5));
    res.push_back({i, S, eS, b.Integral(mD0 - 2 * sig, mD0 + 2 * sig) / bw, h->Integral()});
    h->GetListOfFunctions()->Add(f.Clone("fit"));
  }
  std::sort(res.begin(), res.end(), [](const R &a, const R &b) { return a.S / a.eS > b.S / b.eS; });
  FILE *o = fopen(Form("d0_grid_%.0fMeV.txt", sig * 1000), "w");
  int npos2 = 0, nneg2 = 0;
  for (auto &r : res)
  {
    fprintf(o, "%-45s N=%8.0f S=%7.1f +- %6.1f (%5.2f sigma) B(2s)=%8.0f S/B=%.3f\n", cfgs[r.i].name.c_str(), r.N, r.S, r.eS, r.S / r.eS, r.B, r.S / r.B);
    if (r.S / r.eS > 2) ++npos2;
    if (r.S / r.eS < -2) ++nneg2;
  }
  fprintf(o, "configurations %zu: S/err > +2: %d, < -2: %d\n", res.size(), npos2, nneg2);
  fclose(o);
  printf("configurations %zu: S/err > +2: %d, < -2: %d\n", res.size(), npos2, nneg2);
  TCanvas c("c", "c", 1600, 800);
  c.Divide(4, 2);
  for (int k = 0; k < 8 && k < (int) res.size(); ++k)
  {
    c.cd(k + 1);
    TH1D *h = hs[res[k].i].GetPtr();
    h->SetMinimum(0);
    h->Draw("e");
    ((TF1 *) h->GetListOfFunctions()->FindObject("fit"))->Draw("same");
    TLatex t;
    t.SetNDC();
    t.SetTextSize(0.05);
    t.DrawLatex(0.35, 0.83, Form("S=%.0f#pm%.0f (%.1f#sigma)", res[k].S, res[k].eS, res[k].S / res[k].eS));
  }
  c.SaveAs(Form("d0_grid_best_%.0fMeV.png", sig * 1000));
}
