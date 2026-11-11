// hump2.C : is the pointing/non-pointing hump made by displaced (non-D0-like) vertices and the stored-mass definition?
// usage (note root): root -l -b -q 'scripts/hump2.C+'  -> logs/hump2.txt, logs/hump2_values.txt, figures/hump2_*.{pdf,png}
#include <Math/Vector4D.h>
#include <ROOT/RDataFrame.hxx>
#include <TCanvas.h>
#include <TF1.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TProfile.h>
#include <TStyle.h>

#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

namespace
{
  const char *kFull = "../ntuples/d0tight_full/D0tight_pT1p5_run*.root";
  const double mD0 = 1.86484;
  FILE *gL, *gV;
  template <typename... A>
  void L(const char *f, A... a)
  {
    printf(f, a...);
    fprintf(gL, f, a...);
  }
  void V(const std::string &k, const std::string &v) { fprintf(gV, "%s %s\n", k.c_str(), v.c_str()); }
  TH1D *hratio(const TH1D *a, const TH1D *b, const char *name)
  {
    auto *r = (TH1D *) a->Clone(name);
    r->Reset();
    for (int i = 1; i <= a->GetNbinsX(); ++i)
    {
      double x = a->GetBinContent(i), y = b->GetBinContent(i);
      if (x <= 0 || y <= 0) continue;
      r->SetBinContent(i, x / y);
      r->SetBinError(i, x / y * std::sqrt(1 / x + 1 / y));
    }
    r->Scale(r->GetNbinsX() / r->Integral());
    return r;
  }
  struct Fit
  {
    double dchi, mean, sig, zn, chiP2;
    int ndf;
  };
  Fit fitR(TH1D *r)
  {
    TF1 p2("p2", "pol2", 1.65, 2.08);
    r->Fit(&p2, "QR0");
    TF1 b("b", "[0]*exp(-0.5*((x-[1])/[2])^2)+pol2(3)", 1.65, 2.08);
    b.SetParameters(0.05, 1.88, 0.08, p2.GetParameter(0), p2.GetParameter(1), p2.GetParameter(2));
    b.SetParLimits(1, 1.75, 2.0);
    b.SetParLimits(2, 0.005, 0.25);
    r->Fit(&b, "QR0");
    TF1 n("n", "[0]*exp(-0.5*((x-[1])/[2])^2)+pol2(3)", 1.65, 2.08);
    n.SetParameters(0, mD0, 0.015, p2.GetParameter(0), p2.GetParameter(1), p2.GetParameter(2));
    n.FixParameter(1, mD0);
    n.FixParameter(2, 0.015);
    r->Fit(&n, "QR0");
    return {p2.GetChisquare() - b.GetChisquare(), b.GetParameter(1), b.GetParameter(2), n.GetParError(0) > 0 ? n.GetParameter(0) / n.GetParError(0) : 0,
            p2.GetChisquare(), p2.GetNDF()};
  }
}  // namespace

void hump2()
{
  gL = fopen("logs/hump2.txt", "w");
  gV = fopen("logs/hump2_values.txt", "w");
  ROOT::EnableImplicitMT(6);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetTitleSize(0.055, "XYZ");
  gStyle->SetLabelSize(0.05, "XYZ");
  gStyle->SetPadLeftMargin(0.17);
  gStyle->SetPadBottomMargin(0.14);
  gStyle->SetPadRightMargin(0.04);
  gStyle->SetPadTopMargin(0.06);
  using V4 = ROOT::Math::PtEtaPhiMVector;
  const double mpi = 0.13957, mK = 0.493677, mp = 0.938272, me = 0.000511;
  auto inv = [](double m1, double m2) {
    return [=](float a, float b, float c, float x, float y, float z) { return (float) (V4(a, b, c, m1) + V4(x, y, z, m2)).M(); };
  };
  std::vector<std::string> trk = {"track_1_pT", "track_1_pseudorapidity", "track_1_phi", "track_2_pT", "track_2_pseudorapidity", "track_2_phi"};
  ROOT::RDataFrame d0("DecayTree", kFull);
  auto d = d0.Define("ls", "D0_decayLength/D0_decayLengthErr")
               .Define("mre", inv(mpi, mK), trk)
               .Define("dm", "mre - D0_mass")
               .Define("mpipi", inv(mpi, mpi), trk)
               .Define("mppi", inv(mp, mpi), trk)  // track 1 = proton
               .Define("mpip", inv(mpi, mp), trk)  // track 2 = proton
               .Define("mee", inv(me, me), trk)
               .Define("rxy", "sqrt(D0_x*D0_x+D0_y*D0_y)")
               .Define("Lcm", "D0_decayLength");
  const std::string SIG = "D0_DIRA>0.995 && ls>2", CTL = "D0_DIRA>0.95 && D0_DIRA<0.98";
  auto ds = d.Filter(SIG), dc = d.Filter(CTL);
  auto nS = ds.Count(), nC = dc.Count();

  // ---- composition of the pointing and control samples
  auto fr = [&](ROOT::RDF::RNode n, const char *c) { return n.Filter(c).Count(); };
  std::vector<std::pair<std::string, std::string>> cls = {
      {"Lshort", "Lcm<0.1"},
      {"Lmid", "Lcm>=0.1 && Lcm<1"},
      {"Llong", "Lcm>=1"},
      {"dmBig", "abs(dm)>0.01"},
      {"Ks", "abs(mpipi-0.4976)<0.015"},
      {"Lam", "abs(mppi-1.1157)<0.006 || abs(mpip-1.1157)<0.006"},
      {"Gam", "mee<0.05"},
  };
  std::vector<ROOT::RDF::RResultPtr<ULong64_t>> cS, cC;
  for (auto &c : cls)
  {
    cS.push_back(fr(ds, c.second.c_str()));
    cC.push_back(fr(dc, c.second.c_str()));
  }
  // distributions
  auto hLS = ds.Define("lL", "log10(Lcm)").Histo1D({"hLS", "", 70, -3, 2}, "lL");
  auto hLC = dc.Define("lL", "log10(Lcm)").Histo1D({"hLC", "", 70, -3, 2}, "lL");
  auto hRS = ds.Filter("Lcm>0.1").Histo1D({"hRS", "", 80, 0, 8}, "rxy");
  auto hpp = ds.Histo1D({"hpp", "", 80, 0.3, 0.7}, "mpipi");
  auto hppC = dc.Histo1D({"hppC", "", 80, 0.3, 0.7}, "mpipi");
  auto dmL = d.Define("lL", "log10(Lcm)").Define("adm", "abs(dm)").Profile1D({"dmL", "", 50, -3, 2}, "lL", "adm");
  auto dmH = d.Histo1D({"dmH", "", 100, -0.2, 0.2}, "dm");

  // ---- mass ratios for the variants
  struct Var
  {
    std::string key, label, cut, mass;
  };
  std::vector<Var> vars = {
      {"All", "all candidates, stored mass", "1", "D0_mass"},
      {"Re", "all candidates, mass from stored daughters", "1", "mre"},
      {"Short", "L < 1 mm, stored mass", "Lcm<0.1", "D0_mass"},
      {"Long", "L > 1 mm, stored mass", "Lcm>=0.1", "D0_mass"},
      {"Agree", "|m_{stored}-m_{daughters}| < 10 MeV", "abs(dm)<0.01", "D0_mass"},
      {"Veto", "V0 / #gamma veto", "abs(mpipi-0.4976)>0.015 && abs(mppi-1.1157)>0.006 && abs(mpip-1.1157)>0.006 && mee>0.05", "D0_mass"},
      {"ShortAgree", "L < 1 mm and |#Deltam| < 10 MeV", "Lcm<0.1 && abs(dm)<0.01", "D0_mass"},
  };
  std::vector<ROOT::RDF::RResultPtr<TH1D>> hs, hc;
  std::vector<ROOT::RDF::RResultPtr<ULong64_t>> ns, nc;
  for (auto &v : vars)
  {
    hs.push_back(ds.Filter(v.cut).Histo1D({("s" + v.key).c_str(), "", 43, 1.6, 2.1 - 0.5 / 50 * 7}, v.mass));
    hc.push_back(dc.Filter(v.cut).Histo1D({("c" + v.key).c_str(), "", 43, 1.6, 2.1 - 0.5 / 50 * 7}, v.mass));
    ns.push_back(ds.Filter(v.cut).Count());
    nc.push_back(dc.Filter(v.cut).Count());
  }

  L("pointing (%s) %llu, control (%s) %llu\n", SIG.c_str(), *nS, CTL.c_str(), *nC);
  L("composition (fraction of pointing / control):\n");
  for (size_t i = 0; i < cls.size(); ++i)
  {
    double fs = 100. * *cS[i] / *nS, fc = 100. * *cC[i] / *nC;
    L("  %-6s %-70s %6.2f%% / %6.2f%%\n", cls[i].first.c_str(), cls[i].second.c_str(), fs, fc);
    V("hcomp" + cls[i].first + "S", Form("%.1f", fs));
    V("hcomp" + cls[i].first + "C", Form("%.1f", fc));
  }
  L("|m_stored - m_daughters|: rms %.1f MeV; mean |dm| vs log10(L/cm):\n", 1000 * dmH->GetRMS());
  for (int i = 3; i <= dmL->GetNbinsX(); i += 6) L("  L %.4f cm: <|dm|> %.1f MeV\n", std::pow(10, dmL->GetBinCenter(i)), 1000 * dmL->GetBinContent(i));
  V("hdmShort", Form("%.1f", 1000 * dmL->GetBinContent(dmL->FindBin(-1.5))));
  V("hdmLong", Form("%.0f", 1000 * dmL->GetBinContent(dmL->FindBin(0.5))));
  L("median decay length: pointing %.3f cm, control %.3f cm\n", [&] { double q, p = 0.5; hLS->GetQuantiles(1, &q, &p); return std::pow(10, q); }(),
    [&] { double q, p = 0.5; hLC->GetQuantiles(1, &q, &p); return std::pow(10, q); }());

  L("ratio fits (pol2 vs pol2 + broad gaussian; narrow 15 MeV gaussian at the D0 mass):\n");
  std::vector<TH1D *> rs;
  for (size_t i = 0; i < vars.size(); ++i)
  {
    TH1D *r = hratio(hs[i].GetPtr(), hc[i].GetPtr(), ("r" + vars[i].key).c_str());
    Fit f = fitR(r);
    L("  %-10s %-48s Nsig %8llu Nctl %8llu | pol2 chi2/ndf %6.1f/%d | broad dchi2 %6.1f (mean %.3f sigma %3.0f MeV) | narrow z %5.2f\n", vars[i].key.c_str(),
      vars[i].label.c_str(), *ns[i], *nc[i], f.chiP2, f.ndf, f.dchi, f.mean, 1000 * f.sig, f.zn);
    V("hv" + vars[i].key + "Ns", Form("%llu", *ns[i]));
    V("hv" + vars[i].key + "Chi", Form("%.0f/%d", f.chiP2, f.ndf));
    V("hv" + vars[i].key + "Dchi", Form("%.0f", f.dchi));
    V("hv" + vars[i].key + "Z", Form("%.1f", f.zn));
    rs.push_back(r);
  }

  // ---- figures
  {
    TCanvas c("c", "c", 1500, 520);
    c.Divide(3, 1);
    c.cd(1);
    gPad->SetLogy();
    hLS->Scale(1. / hLS->Integral());
    hLC->Scale(1. / hLC->Integral());
    hLS->SetLineColor(kRed + 1);
    hLS->SetLineWidth(2);
    hLC->SetLineWidth(2);
    hLC->GetXaxis()->SetTitle("log_{10}(decay length / cm)");
    hLC->GetYaxis()->SetTitle("fraction");
    hLC->SetMaximum(0.2);
    hLC->DrawCopy("hist");
    hLS->DrawCopy("hist same");
    TLine l(-1, 1e-5, -1, 0.2);
    l.SetLineStyle(2);
    l.DrawClone();
    TLegend lg(0.2, 0.78, 0.7, 0.92);
    lg.SetBorderSize(0);
    lg.SetTextSize(0.045);
    lg.AddEntry(hLS.GetPtr(), "pointing", "l");
    lg.AddEntry(hLC.GetPtr(), "non-pointing", "l");
    lg.DrawClone();
    c.cd(2);
    dmL->GetXaxis()->SetTitle("log_{10}(decay length / cm)");
    dmL->GetYaxis()->SetTitle("#LT|m_{stored} - m_{daughters}|#GT (GeV/c^{2})");
    dmL->SetMarkerStyle(20);
    dmL->SetMarkerSize(0.6);
    dmL->DrawCopy("e");
    c.cd(3);
    hRS->GetXaxis()->SetTitle("vertex radius r_{xy} (cm), pointing, L > 1 mm");
    hRS->GetYaxis()->SetTitle("candidates");
    hRS->SetLineWidth(2);
    hRS->DrawCopy("hist");
    c.SaveAs("figures/hump2_displaced.pdf");
    c.SaveAs("figures/hump2_displaced.png");
  }
  {
    TCanvas c("c", "c", 1500, 520);
    c.Divide(3, 1);
    int col[7] = {kBlack, kGray + 2, kBlue + 1, kRed + 1, kGreen + 2, kMagenta + 1, kOrange + 7};
    std::vector<std::vector<int>> panels = {{0, 2, 3}, {0, 1, 4}, {0, 5, 6}};
    for (int p = 0; p < 3; ++p)
    {
      c.cd(p + 1);
      TLegend lg(0.2, 0.74, 0.95, 0.93);
      lg.SetBorderSize(0);
      lg.SetTextSize(0.04);
      for (size_t k = 0; k < panels[p].size(); ++k)
      {
        int i = panels[p][k];
        TH1D *r = rs[i];
        r->SetMarkerStyle(20 + k);
        r->SetMarkerSize(0.7);
        r->SetMarkerColor(col[i]);
        r->SetLineColor(col[i]);
        r->GetXaxis()->SetTitle("m(K#pi) (GeV/c^{2})");
        r->GetYaxis()->SetTitle("pointing / non-pointing (norm.)");
        r->SetMinimum(0.6);
        r->SetMaximum(1.45);
        r->DrawCopy(k ? "e same" : "e");
        lg.AddEntry(r, vars[i].label.c_str(), "pe");
      }
      lg.DrawClone();
      TLine l(mD0, 0.6, mD0, 1.45);
      l.SetLineStyle(3);
      l.DrawClone();
    }
    c.SaveAs("figures/hump2_ratios.pdf");
    c.SaveAs("figures/hump2_ratios.png");
  }
  fclose(gL);
  fclose(gV);
}
