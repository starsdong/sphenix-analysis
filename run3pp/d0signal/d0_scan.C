// d0_scan.C : scan topological / PID selections on the D0 skim and fit for a D0 -> K pi peak
// usage: root -b -q d0_scan.C+            output: d0_scan.txt, d0_scan_*.png
#include <Math/Vector4D.h>
#include <Math/Boost.h>
#include <ROOT/RDataFrame.hxx>
#include <TCanvas.h>
#include <TF1.h>
#include <TGraph.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLatex.h>
#include <TStyle.h>

#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

namespace
{
  // pion dE/dx band: median and width vs p, measured from the pion-hypothesis daughter (pion dominated)
  TGraph gMed, gWid;
  double piMed(double p) { return gMed.Eval(std::min(std::max(p, 0.55), 2.95)); }
  double piWid(double p) { return gWid.Eval(std::min(std::max(p, 0.55), 2.95)); }

  struct FitRes
  {
    double S, eS, B, mean, sigma, chi2ndf;
  };
  FitRes fitPeak(TH1D *h, bool fixSigma = false, double sig0 = 0.015)
  {
    FitRes r{0, 0, 0, 0, 0, 0};
    if (h->GetEntries() < 50) return r;
    TF1 f("f", "gaus(0)+pol2(3)", 1.70, 2.05);
    // background-only start
    TF1 fb("fb", "pol2", 1.70, 2.05);
    h->Fit(&fb, "QN0R");
    f.SetParameters(h->GetMaximum() * 0.05, 1.865, sig0, fb.GetParameter(0), fb.GetParameter(1), fb.GetParameter(2));
    f.SetParLimits(0, 0, h->GetMaximum());
    f.SetParLimits(1, 1.845, 1.885);
    if (fixSigma)
      f.FixParameter(2, sig0);
    else
      f.SetParLimits(2, 0.006, 0.035);
    h->Fit(&f, "QN0RL");
    double bw = h->GetBinWidth(1);
    r.S = f.GetParameter(0) * f.GetParameter(2) * std::sqrt(2 * M_PI) / bw;
    double rel = f.GetParameter(0) > 0 ? f.GetParError(0) / f.GetParameter(0) : 1;
    r.eS = r.S * rel;
    r.mean = f.GetParameter(1);
    r.sigma = f.GetParameter(2);
    TF1 bonly("bonly", "pol2", 1.70, 2.05);
    bonly.SetParameters(f.GetParameter(3), f.GetParameter(4), f.GetParameter(5));
    r.B = bonly.Integral(r.mean - 2 * r.sigma, r.mean + 2 * r.sigma) / bw;
    r.chi2ndf = f.GetNDF() ? f.GetChisquare() / f.GetNDF() : 0;
    h->GetListOfFunctions()->Add(f.Clone("fit"));
    return r;
  }
}  // namespace

void d0_scan()
{
  ROOT::EnableImplicitMT(6);
  gStyle->SetOptStat(0);
  ROOT::RDataFrame d0("DecayTree", "../ntuples/d0skim/D0skim_prod10pct_*.root");

  // ---- 1) pion dE/dx band from the pion-hypothesis daughter
  auto d1 = d0.Define("pPi", "abs(track_1_PDG_ID)==211 ? track_1_pT*cosh(track_1_pseudorapidity) : track_2_pT*cosh(track_2_pseudorapidity)")
                .Define("dPi", "abs(track_1_PDG_ID)==211 ? track_1_calculated_dEdx : track_2_calculated_dEdx");
  auto hband = d1.Histo2D({"hband", "", 49, 0.5, 2.95, 600, 0, 1500}, "pPi", "dPi");
  auto nK2 = d0.Filter("abs(track_2_PDG_ID)==321").Count();
  auto nAll = d0.Count();
  int ip = 0;
  for (int b = 1; b <= hband->GetNbinsX(); ++b)
  {
    auto *py = hband->ProjectionY("py", b, b);
    if (py->GetEntries() < 1000) continue;
    double q[3], pr[3] = {0.25, 0.5, 0.75};
    py->GetQuantiles(3, q, pr);
    // width from the lower half (the upper half is contaminated by K / p)
    gMed.SetPoint(ip, hband->GetXaxis()->GetBinCenter(b), q[1]);
    gWid.SetPoint(ip, hband->GetXaxis()->GetBinCenter(b), (q[1] - q[0]) / 0.674);
    ++ip;
  }
  printf("candidates %llu, track_2 is the kaon in %.1f%%\n", *nAll, 100. * *nK2 / *nAll);
  for (double p : {0.55, 0.7, 0.9, 1.2, 1.6, 2.2}) printf("pion band p=%.2f: median %.0f width %.0f\n", p, piMed(p), piWid(p));

  // ---- 2) derived variables
  auto d = d0.Define("kIs1", "abs(track_1_PDG_ID)==321")
               .Define("pK", "kIs1 ? track_1_pT*cosh(track_1_pseudorapidity) : track_2_pT*cosh(track_2_pseudorapidity)")
               .Define("pP", "kIs1 ? track_2_pT*cosh(track_2_pseudorapidity) : track_1_pT*cosh(track_1_pseudorapidity)")
               .Define("dK", "kIs1 ? track_1_calculated_dEdx : track_2_calculated_dEdx")
               .Define("dP", "kIs1 ? track_2_calculated_dEdx : track_1_calculated_dEdx")
               .Define("zK", [](float p, float de) { return (de - piMed(p)) / piWid(p); }, {"pK", "dK"})
               .Define("zP", [](float p, float de) { return (de - piMed(p)) / piWid(p); }, {"pP", "dP"})
               .Define("ls", "D0_decayLength/D0_decayLengthErr")
               .Define("lsxy", "D0_decayLength_xy/D0_decayLengthErr_xy")
               .Define("minDCAsig", "std::min(abs(track_1_PV_DCA_sig),abs(track_2_PV_DCA_sig))")
               .Define("prodDCA", "track_1_PV_DCA_xy*track_2_PV_DCA_xy")
               .Define("cts",
                       [](bool kIs1, float pt1, float eta1, float phi1, float pt2, float eta2, float phi2) {
                         using V = ROOT::Math::PtEtaPhiMVector;
                         V k = kIs1 ? V(pt1, eta1, phi1, 0.493677) : V(pt2, eta2, phi2, 0.493677);
                         V pi = kIs1 ? V(pt2, eta2, phi2, 0.139570) : V(pt1, eta1, phi1, 0.139570);
                         auto D = k + pi;
                         ROOT::Math::Boost bst(D.BoostToCM());
                         auto kr = bst(k);
                         return (double) (kr.Vect().Unit().Dot(D.Vect().Unit()));
                       },
                       {"kIs1", "track_1_pT", "track_1_pseudorapidity", "track_1_phi", "track_2_pT", "track_2_pseudorapidity", "track_2_phi"})
               // kaon PID: below 0.9 GeV/c the kaon must lie above the pion band; pion must be pion-like (or high p)
               .Define("kPID", "pK > 0.9 || zK > 2.0")
               .Define("kPIDstrict", "(pK < 0.9 && zK > 2.0) || (pK >= 0.9 && pK < 1.2 && zK > 1.0)")
               .Define("piPID", "abs(zP) < 3");

  std::string topoA = "D0_DIRA>0.99 && ls>3";
  std::string topoB = "D0_DIRA>0.995 && ls>4 && track_1_track_2_DCA<0.01";
  std::string topoC = "D0_DIRA>0.999 && ls>5 && track_1_track_2_DCA<0.008 && minDCAsig>2.5";
  std::vector<std::pair<std::string, std::string>> sel = {
      {"00 skim (all)", "1"},
      {"01 DIRA>0.99, L/sL>3", topoA},
      {"02 topoB: DIRA>0.995, L/sL>4, DCA12<0.01", topoB},
      {"03 topoC: DIRA>0.999, L/sL>5, DCA12<80um, minDCAsig>2.5", topoC},
      {"04 topoB + |cts|<0.8", topoB + " && abs(cts)<0.8"},
      {"05 topoB + PID", topoB + " && kPID && piPID"},
      {"06 topoB + cts + PID", topoB + " && abs(cts)<0.8 && kPID && piPID"},
      {"07 topoB + cts + PID, pT>1", topoB + " && abs(cts)<0.8 && kPID && piPID && D0_pT>1"},
      {"08 topoB + cts + PID, pT>2", topoB + " && abs(cts)<0.8 && kPID && piPID && D0_pT>2"},
      {"09 topoC + cts + PID", topoC + " && abs(cts)<0.8 && kPID && piPID"},
      {"10 topoC + cts + PID, pT>1.5", topoC + " && abs(cts)<0.8 && kPID && piPID && D0_pT>1.5"},
      {"11 topoB + cts + strict K PID (pK<1.2)", topoB + " && abs(cts)<0.8 && kPIDstrict && piPID"},
      {"12 topoA + cts + strict K PID (pK<1.2)", topoA + " && abs(cts)<0.8 && kPIDstrict && piPID"},
      {"13 topoC + strict K PID", topoC + " && kPIDstrict && piPID"},
      {"14 topoB + D0 PV DCA<50um + cts + PID", topoB + " && D0_PV_DCA<0.005 && abs(cts)<0.8 && kPID && piPID"},
      {"15 topoB + prodDCA<0 + cts + PID", topoB + " && prodDCA<0 && abs(cts)<0.8 && kPID && piPID"},
  };
  std::vector<ROOT::RDF::RResultPtr<TH1D>> hs;
  std::vector<ROOT::RDF::RResultPtr<ULong64_t>> ns;
  for (size_t i = 0; i < sel.size(); ++i)
  {
    auto f = d.Filter(sel[i].second);
    hs.push_back(f.Histo1D({Form("hm%zu", i), (sel[i].first + ";m(K#pi) (GeV/c^{2});counts / 5 MeV").c_str(), 100, 1.6, 2.1}, "D0_mass"));
    ns.push_back(f.Count());
  }
  // control distributions for PID
  auto hzK = d.Filter("pK<0.9").Histo2D({"hzK", "kaon-hyp daughter, p<0.9;p (GeV/c);z (pion band)", 40, 0.5, 0.9, 120, -4, 20}, "pK", "zK");

  FILE *out = fopen("d0_scan.txt", "w");
  TCanvas c("c", "c", 1800, 1600);
  c.Divide(4, 4);
  for (size_t i = 0; i < sel.size(); ++i)
  {
    TH1D *h = hs[i].GetPtr();
    FitRes r = fitPeak(h);
    FitRes rf = fitPeak((TH1D *) h->Clone(Form("hc%zu", i)), true, 0.015);
    double signif = r.S > 0 ? r.S / std::sqrt(r.S + r.B) : 0;
    char line[512];
    snprintf(line, sizeof line,
             "%-58s N=%9llu | free: S=%7.0f+-%5.0f mean=%.4f sigma=%.4f S/sqrt(S+B)=%5.2f S/B=%.3f chi2/ndf=%.2f | sigma fixed 15 MeV: S=%7.0f+-%5.0f mean=%.4f\n",
             sel[i].first.c_str(), *ns[i], r.S, r.eS, r.mean, r.sigma, signif, r.B > 0 ? r.S / r.B : 0, r.chi2ndf, rf.S, rf.eS, rf.mean);
    printf("%s", line);
    fputs(line, out);
    c.cd(i + 1);
    h->SetMinimum(0);
    h->Draw("e");
    if (auto *fit = h->GetListOfFunctions()->FindObject("fit")) ((TF1 *) fit)->Draw("same");
    TLatex t;
    t.SetNDC();
    t.SetTextSize(0.05);
    t.DrawLatex(0.4, 0.83, Form("S=%.0f#pm%.0f  #sigma=%.0f MeV", r.S, r.eS, 1000 * r.sigma));
  }
  fclose(out);
  c.SaveAs("d0_scan_mass.png");
  TCanvas c2("c2", "c2", 700, 500);
  c2.SetLogz();
  hzK->Draw("colz");
  c2.SaveAs("d0_scan_zK.png");
}
