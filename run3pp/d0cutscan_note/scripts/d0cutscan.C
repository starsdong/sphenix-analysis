// d0cutscan.C : all numbers, tables and figures of the D0 cut-scan note.
// usage (from the note root):  root -l -b -q 'scripts/d0cutscan.C+("<step>")'   step = inputs | scan | grid | trial | ratio
// Inputs (local copies, see logs/local_files_md5.txt):
//   10% skim : ../ntuples/d0skim/D0skim_prod10pct_*.root
//   full skim with D0 pT>1.5 : ../ntuples/d0tight_full/D0tight_pT1p5_run*.root
// Outputs: logs/<step>.txt (human readable), logs/<step>_values.txt ("KEY value", turned into LaTeX macros by
// scripts/make_numbers.py), tables/<step>_*.tex, figures/<name>.{pdf,png}
#include <Math/Boost.h>
#include <Math/Vector4D.h>
#include <ROOT/RDataFrame.hxx>
#include <TCanvas.h>
#include <TF1.h>
#include <TGraph.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TPad.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TSystem.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

namespace
{
  const char *k10 = "../ntuples/d0skim/D0skim_prod10pct_*.root";
  const char *kFull = "../ntuples/d0tight_full/D0tight_pT1p5_run*.root";
  const double mD0 = 1.86484;  // PDG
  const double kLo = 1.60, kHi = 2.10, kFitLo = 1.70, kFitHi = 2.05;
  const int kNb = 100;  // 5 MeV bins

  FILE *gLog = nullptr, *gVal = nullptr;
  void openLogs(const std::string &step)
  {
    gLog = fopen(("logs/" + step + ".txt").c_str(), "w");
    gVal = fopen(("logs/" + step + "_values.txt").c_str(), "w");
  }
  void closeLogs()
  {
    fclose(gLog);
    fclose(gVal);
  }
  template <typename... A>
  void L(const char *fmt, A... a)
  {
    printf(fmt, a...);
    fprintf(gLog, fmt, a...);
  }
  void V(const std::string &k, const std::string &v) { fprintf(gVal, "%s %s\n", k.c_str(), v.c_str()); }
  std::string f0(double x) { return Form("%.0f", x); }
  std::string f1(double x) { return Form("%.1f", x); }
  std::string f2(double x) { return Form("%.2f", x); }
  std::string f3(double x) { return Form("%.3f", x); }
  std::string big(double x)  // 1234567 -> 1{,}234{,}567
  {
    std::string s = Form("%.0f", x), o;
    int n = s.size();
    for (int i = 0; i < n; ++i)
    {
      o += s[i];
      if ((n - 1 - i) % 3 == 0 && i != n - 1) o += "{,}";
    }
    return o;
  }

  void style()
  {
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetTitleSize(0.055, "XYZ");
    gStyle->SetLabelSize(0.05, "XYZ");
    gStyle->SetTitleOffset(1.1, "Y");
    gStyle->SetPadLeftMargin(0.15);
    gStyle->SetPadBottomMargin(0.14);
    gStyle->SetPadRightMargin(0.04);
    gStyle->SetPadTopMargin(0.06);
  }
  void save(TCanvas &c, const std::string &name)
  {
    c.SaveAs(("figures/" + name + ".pdf").c_str());
    c.SaveAs(("figures/" + name + ".png").c_str());
  }

  // ---------------- pion dE/dx band (median, width from the 25%..50% quantiles) of the pion-hypothesis daughter
  TGraph gMed, gWid;
  void loadBand()
  {
    std::ifstream in("logs/pionband.txt");
    std::string line;
    int i = 0;
    while (std::getline(in, line))
    {
      if (line.empty() || line[0] == '#') continue;
      double p, m, w;
      if (sscanf(line.c_str(), "%lf %lf %lf", &p, &m, &w) == 3)
      {
        gMed.SetPoint(i, p, m);
        gWid.SetPoint(i, p, w);
        ++i;
      }
    }
    if (i < 3)
    {
      printf("run step inputs first (logs/pionband.txt missing)\n");
      gSystem->Exit(1);
    }
  }
  double pmin() { return gMed.GetX()[0]; }
  double pmax() { return gMed.GetX()[gMed.GetN() - 1]; }
  float zPi(float p, float de)
  {
    double pp = std::min(std::max((double) p, pmin()), pmax());
    return (de - gMed.Eval(pp)) / gWid.Eval(pp);
  }

  // common derived columns (track_1 = pion hypothesis, track_2 = kaon hypothesis; checked in step inputs)
  ROOT::RDF::RNode derive(ROOT::RDF::RNode d)
  {
    return d.Define("pK", "track_2_pT*cosh(track_2_pseudorapidity)")
        .Define("pP", "track_1_pT*cosh(track_1_pseudorapidity)")
        .Define("zK", [](float p, float de) { return zPi(p, de); }, {"pK", "track_2_calculated_dEdx"})
        .Define("zP", [](float p, float de) { return zPi(p, de); }, {"pP", "track_1_calculated_dEdx"})
        .Define("ls", "D0_decayLength/D0_decayLengthErr")
        .Define("minDCAsig", "std::min(abs(track_1_PV_DCA_sig),abs(track_2_PV_DCA_sig))")
        .Define("cts",
                [](float pt1, float eta1, float phi1, float pt2, float eta2, float phi2) {
                  using V4 = ROOT::Math::PtEtaPhiMVector;
                  V4 pi(pt1, eta1, phi1, 0.139570), k(pt2, eta2, phi2, 0.493677);
                  auto D = k + pi;
                  auto kr = ROOT::Math::Boost(D.BoostToCM())(k);
                  return (double) kr.Vect().Unit().Dot(D.Vect().Unit());
                },
                {"track_1_pT", "track_1_pseudorapidity", "track_1_phi", "track_2_pT", "track_2_pseudorapidity", "track_2_phi"})
        .Define("kPID", "pK > 0.9 || zK > 2.0")
        .Define("piPID", "abs(zP) < 3")
        .Define("kPIDstrict", "(pK < 0.9 && zK > 2.0) || (pK >= 0.9 && pK < 1.2 && zK > 1.0)");
  }

  // ---------------- fits on a mass histogram (counts): gaus at fixed mean/sigma + pol2, amplitude unbounded
  struct Fix
  {
    double S, eS, B, chi2, ndf;
  };
  Fix fitFixed(TH1D *h, double mean, double sig, bool keep = false, const char *opt = "QN0RL")
  {
    TF1 f(Form("ffix_%s", h->GetName()), "[0]*exp(-0.5*((x-[1])/[2])^2)+pol2(3)", kFitLo, kFitHi);
    TF1 b0("b0", "pol2", kFitLo, kFitHi);
    h->Fit(&b0, opt);
    f.SetParameters(0, mean, sig, b0.GetParameter(0), b0.GetParameter(1), b0.GetParameter(2));
    f.FixParameter(1, mean);
    f.FixParameter(2, sig);
    h->Fit(&f, opt);
    double bw = h->GetBinWidth(1), k = sig * std::sqrt(2 * M_PI) / bw;
    TF1 b("b", "pol2", kFitLo, kFitHi);
    b.SetParameters(f.GetParameter(3), f.GetParameter(4), f.GetParameter(5));
    Fix r{f.GetParameter(0) * k, f.GetParError(0) * k, b.Integral(mean - 2 * sig, mean + 2 * sig) / bw, f.GetChisquare(), (double) f.GetNDF()};
    if (keep) h->GetListOfFunctions()->Add(f.Clone("fitfix"));
    return r;
  }
  struct Free
  {
    double S, eS, mean, sig;
    bool atLimit;
  };
  Free fitFree(TH1D *h)
  {
    TF1 f(Form("ffree_%s", h->GetName()), "[0]*exp(-0.5*((x-[1])/[2])^2)+pol2(3)", kFitLo, kFitHi);
    TF1 b0("b0", "pol2", kFitLo, kFitHi);
    h->Fit(&b0, "QN0RL");
    f.SetParameters(h->GetMaximum() * 0.02, mD0, 0.015, b0.GetParameter(0), b0.GetParameter(1), b0.GetParameter(2));
    f.SetParLimits(0, 0, h->GetMaximum());
    f.SetParLimits(1, 1.845, 1.885);
    f.SetParLimits(2, 0.006, 0.035);
    h->Fit(&f, "QN0RL");
    double bw = h->GetBinWidth(1), k = f.GetParameter(2) * std::sqrt(2 * M_PI) / bw;
    auto near = [](double x, double lo, double hi) { return std::fabs(x - lo) < 1e-3 * (hi - lo) || std::fabs(x - hi) < 1e-3 * (hi - lo); };
    bool lim = near(f.GetParameter(1), 1.845, 1.885) || near(f.GetParameter(2), 0.006, 0.035) || f.GetParameter(0) < 1e-6 * h->GetMaximum();
    return {f.GetParameter(0) * k, f.GetParError(0) * k, f.GetParameter(1), f.GetParameter(2), lim};
  }
  // sideband cross-check: pol2 fitted to the fit range with the window mean+-3 sigma excluded; signal = counts in
  // mean+-2 sigma minus the pol2 integral there; error = sqrt(counts in the window) (background-fit error neglected)
  std::pair<double, double> sideband(TH1D *h, double mean, double sig)
  {
    TF1 f(Form("fsb_%s", h->GetName()),
          [=](double *x, double *p) {
            if (std::fabs(x[0] - mean) < 3 * sig) TF1::RejectPoint();
            return p[0] + p[1] * x[0] + p[2] * x[0] * x[0];
          },
          kFitLo, kFitHi, 3);
    h->Fit(&f, "QN0RL");
    int b1 = h->FindBin(mean - 2 * sig + 1e-6), b2 = h->FindBin(mean + 2 * sig - 1e-6);
    double nw = h->Integral(b1, b2), bw = h->GetBinWidth(1);
    double lo = h->GetBinLowEdge(b1), hi = h->GetBinLowEdge(b2) + bw;
    double nb = (f.GetParameter(0) * (hi - lo) + f.GetParameter(1) * (hi * hi - lo * lo) / 2 + f.GetParameter(2) * (hi * hi * hi - lo * lo * lo) / 3) / bw;
    return {nw - nb, std::sqrt(nw)};
  }
  void drawMass(TH1D *h, const std::string &label, double S, double eS, bool fit = true)
  {
    h->SetMarkerStyle(20);
    h->SetMarkerSize(0.5);
    h->SetLineColor(kBlack);
    h->SetMinimum(0);
    h->SetMaximum(h->GetMaximum() * 1.35);
    h->GetXaxis()->SetTitle("m(K#pi) (GeV/c^{2})");
    h->GetXaxis()->SetNdivisions(505);
    h->GetYaxis()->SetTitle("candidates / 5 MeV");
    h->GetYaxis()->SetMaxDigits(3);
    h->GetYaxis()->SetTitleOffset(1.35);
    gPad->SetLeftMargin(0.17);
    h->Draw("e");
    if (fit)
      if (auto *f = (TF1 *) h->GetListOfFunctions()->FindObject("fitfix"))
      {
        f->SetLineColor(kRed + 1);
        f->SetLineWidth(2);
        f->Draw("same");
      }
    TLine l(mD0, 0, mD0, h->GetMaximum());
    l.SetLineStyle(3);
    l.SetLineColor(kGray + 2);
    l.DrawClone();
    TLatex t;
    t.SetNDC();
    t.SetTextSize(0.052);
    t.DrawLatex(0.18, 0.86, label.c_str());
    if (fit) t.DrawLatex(0.18, 0.79, Form("S(m_{D^{0}}, #sigma=15 MeV) = %.0f #pm %.0f", S, eS));
  }
}  // namespace

// =====================================================================================================
// step inputs: sample content, mass-hypothesis assignment, duplicate pairs, pion dE/dx band, dE/dx figure
void stepInputs()
{
  openLogs("inputs");
  for (auto [tag, files] : std::vector<std::pair<std::string, const char *>>{{"Ten", k10}, {"Full", kFull}})
  {
    ROOT::RDataFrame d("DecayTree", files);
    auto n = d.Count();
    auto nK2 = d.Filter("abs(track_2_PDG_ID)==321 && abs(track_1_PDG_ID)==211").Count();
    auto nOpp = d.Filter("int(track_1_charge)*int(track_2_charge)<0").Count();
    auto nD0 = d.Filter("D0_PDG_ID==421").Count();
    auto nKminus = d.Filter("(D0_PDG_ID==421 && track_2_charge<0) || (D0_PDG_ID==-421 && track_2_charge>0)").Count();
    auto mmin = d.Min("D0_mass"), mmax = d.Max("D0_mass");
    auto ptmin = d.Min("D0_pT");
    auto diramin = d.Min("D0_DIRA");
    auto t1pt = d.Min("track_1_pT"), t2pt = d.Min("track_2_pT");
    auto dca1 = d.Define("a", "abs(track_1_PV_DCA)").Min("a"), dca2 = d.Define("a", "abs(track_2_PV_DCA)").Min("a");
    L("[%s] %s\n  candidates %llu\n  track_1 = pi hypothesis and track_2 = K hypothesis: %llu (%.4f%%)\n", tag.c_str(), files, *n, *nK2, 100. * *nK2 / *n);
    L("  opposite-sign daughters: %llu\n  D0 (vs D0bar) %llu; kaon charge consistent with D0/D0bar: %llu\n", *nOpp, *nD0, *nKminus);
    L("  ranges: mass %.4f-%.4f, D0 pT min %.3f, DIRA min %.4f, daughter pT min %.3f %.3f, |PV DCA| min %.4f %.4f cm\n", *mmin, *mmax, *ptmin, *diramin,
      *t1pt, *t2pt, *dca1, *dca2);
    V("n" + tag, big(*n));
    V("hypFrac" + tag, Form("%.2f", 100. * *nK2 / *n));
    V("ptMin" + tag, f3(*ptmin));
    V("diraMin" + tag, Form("%.4f", *diramin));

    // duplicate track pairs: same run/event, same two track IDs (any role)
    auto run = d.Take<int>("runNumber");
    auto ev = d.Take<int>("eventNumber");
    auto id1 = d.Take<int>("track_1_track_ID");
    auto id2 = d.Take<int>("track_2_track_ID");
    // sort (run, event, min id, max id, id of the pion-hypothesis track) and walk the groups
    std::vector<std::array<int, 5>> k(run->size());
    for (size_t i = 0; i < run->size(); ++i)
    {
      int a = (*id1)[i], b = (*id2)[i];
      k[i] = {(*run)[i], (*ev)[i], std::min(a, b), std::max(a, b), a};
    }
    std::sort(k.begin(), k.end());
    long nPairs = 0, nMulti = 0, nSwap = 0, nSamePVdiff = 0;
    for (size_t i = 0; i < k.size();)
    {
      size_t j = i;
      std::set<int> roles;
      while (j < k.size() && k[j][0] == k[i][0] && k[j][1] == k[i][1] && k[j][2] == k[i][2] && k[j][3] == k[i][3]) roles.insert(k[j++][4]);
      ++nPairs;
      if (j - i > 1)
      {
        ++nMulti;
        if (roles.size() > 1) ++nSwap;  // the same pair stored with K and pi exchanged
        else ++nSamePVdiff;             // the same assignment stored more than once (e.g. other PV)
      }
      i = j;
    }
    L("  unique track pairs %ld; stored more than once %ld (%.3f%%): with K<->pi exchanged %ld, same assignment (other PV) %ld\n", nPairs, nMulti,
      100. * nMulti / nPairs, nSwap, nSamePVdiff);
    V("pairs" + tag, big(nPairs));
    V("multi" + tag, big(nMulti));
    V("multiFrac" + tag, Form("%.3f", 100. * nMulti / nPairs));
    V("swap" + tag, big(nSwap));
    V("samePV" + tag, big(nSamePVdiff));
  }

  // pion band from the pion-hypothesis daughter of the 10% sample
  style();
  ROOT::RDataFrame d("DecayTree", k10);
  auto dd = d.Define("pP", "track_1_pT*cosh(track_1_pseudorapidity)").Define("pK", "track_2_pT*cosh(track_2_pseudorapidity)");
  auto hb = dd.Histo2D({"hb", "", 34, 0.5, 2.2, 600, 0, 1500}, "pP", "track_1_calculated_dEdx");
  auto h1 = dd.Histo2D({"h1", ";p (GeV/c);dE/dx (a.u.)", 100, 0.4, 3.0, 150, 0, 1200}, "pP", "track_1_calculated_dEdx");
  auto h2 = dd.Histo2D({"h2", ";p (GeV/c);dE/dx (a.u.)", 100, 0.4, 3.0, 150, 0, 1200}, "pK", "track_2_calculated_dEdx");
  FILE *pb = fopen("logs/pionband.txt", "w");
  fprintf(pb, "# p(GeV/c) median width ; pion-hypothesis daughter, 10%% skim; width = (median - q25)/0.674\n");
  FILE *tb = fopen("tables/pionband.tex", "w");
  L("pion band (p, median, width, entries):\n");
  for (int b = 1; b <= hb->GetNbinsX(); ++b)
  {
    auto *py = hb->ProjectionY("py", b, b);
    if (py->GetEntries() < 2000) continue;
    double q[2], pr[2] = {0.25, 0.5};
    py->GetQuantiles(2, q, pr);
    double p = hb->GetXaxis()->GetBinCenter(b);
    fprintf(pb, "%.3f %.1f %.1f\n", p, q[1], (q[1] - q[0]) / 0.674);
    L("  %.3f %.1f %.1f %.0f\n", p, q[1], (q[1] - q[0]) / 0.674, py->GetEntries());
    if (b % 4 == 1) fprintf(tb, "%.2f & %.0f & %.0f & %s \\\\\n", p, q[1], (q[1] - q[0]) / 0.674, big(py->GetEntries()).c_str());
  }
  fclose(pb);
  fclose(tb);
  loadBand();
  V("bandMedLo", f0(gMed.Eval(0.55)));
  V("bandMedHi", f0(gMed.Eval(2.0)));
  V("bandWidLo", f0(gWid.Eval(0.55)));
  V("bandWidHi", f0(gWid.Eval(2.0)));
  // kaon-hypothesis daughter below 0.9 GeV/c: fraction above the pion band (z>2)
  auto dz = derive(ROOT::RDF::RNode(d));
  auto nlow = dz.Filter("pK<0.9").Count(), nlowK = dz.Filter("pK<0.9 && zK>2").Count(), nAll = dz.Count();
  auto npi = dz.Filter("abs(zP)<3").Count(), nkpid = dz.Filter("kPID && piPID").Count();
  L("kaon-hyp daughter p<0.9: %llu of %llu (%.1f%%); of those z>2: %llu (%.1f%%); pion |z|<3: %.1f%%; kPID&&piPID: %.1f%%\n", *nlow, *nAll,
    100. * *nlow / *nAll, *nlowK, 100. * *nlowK / *nlow, 100. * *npi / *nAll, 100. * *nkpid / *nAll);
  V("kLowFrac", f1(100. * *nlow / *nAll));
  V("kLowAbove", f1(100. * *nlowK / *nlow));
  V("pidKeep", f1(100. * *nkpid / *nAll));

  TCanvas c("c", "c", 1400, 600);
  c.Divide(2, 1);
  for (int i = 0; i < 2; ++i)
  {
    c.cd(i + 1);
    gPad->SetLogz();
    gPad->SetRightMargin(0.12);
    TH2D *h = i == 0 ? h1.GetPtr() : h2.GetPtr();
    h->DrawCopy("colz");
    TGraph up, mid;
    for (int k = 0; k < gMed.GetN(); ++k)
    {
      mid.SetPoint(k, gMed.GetX()[k], gMed.GetY()[k]);
      up.SetPoint(k, gMed.GetX()[k], gMed.GetY()[k] + 2 * gWid.GetY()[k]);
    }
    mid.SetLineColor(kBlack);
    mid.SetLineWidth(2);
    up.SetLineColor(kRed + 1);
    up.SetLineWidth(2);
    up.SetLineStyle(2);
    mid.DrawClone("L");
    up.DrawClone("L");
    TLatex t;
    t.SetNDC();
    t.SetTextSize(0.05);
    t.DrawLatex(0.4, 0.86, i == 0 ? "track 1 (#pi hypothesis)" : "track 2 (K hypothesis)");
    if (i == 1)
    {
      TLegend lg(0.4, 0.68, 0.86, 0.82);
      lg.SetBorderSize(0);
      lg.SetTextSize(0.042);
      lg.AddEntry(&mid, "pion band median", "l");
      lg.AddEntry(&up, "median + 2 width (z = 2)", "l");
      lg.DrawClone();
    }
  }
  save(c, "dedx_bands");
  closeLogs();
}

// ---------------- selections of the hand scan
struct Sel
{
  std::string key, label, cut;
};
std::vector<Sel> scanSelections()
{
  std::string A = "D0_DIRA>0.99 && ls>3";
  std::string B = "D0_DIRA>0.995 && ls>4 && track_1_track_2_DCA<0.01";
  std::string C = "D0_DIRA>0.999 && ls>5 && track_1_track_2_DCA<0.008 && minDCAsig>2.5";
  std::string pc = " && abs(cts)<0.8 && kPID && piPID";
  return {
      {"S0", "skim (no further cut)", "1"},
      {"S1", "A", A},
      {"S2", "B", B},
      {"S3", "C", C},
      {"S4", "B + $|\\cos\\theta^*|<0.8$", B + " && abs(cts)<0.8"},
      {"S5", "B + PID", B + " && kPID && piPID"},
      {"S6", "B + cts + PID", B + pc},
      {"S7", "B + cts + PID, $p_T>1$", B + pc + " && D0_pT>1"},
      {"S8", "B + cts + PID, $p_T>2$", B + pc + " && D0_pT>2"},
      {"S9", "C + cts + PID", C + pc},
      {"S10", "C + cts + PID, $p_T>1.5$", C + pc + " && D0_pT>1.5"},
      {"S11", "B + cts + strict K PID", B + " && abs(cts)<0.8 && kPIDstrict && piPID"},
      {"S12", "A + cts + strict K PID", A + " && abs(cts)<0.8 && kPIDstrict && piPID"},
      {"S13", "C + strict K PID", C + " && kPIDstrict && piPID"},
      {"S14", "B + $D^0$ PV DCA$<50\\,\\mu$m + cts + PID", B + " && D0_PV_DCA<0.005" + pc},
      {"S15", "B + $d_1^{xy} d_2^{xy}<0$ + cts + PID", B + " && track_1_PV_DCA_xy*track_2_PV_DCA_xy<0" + pc},
  };
}

void stepScan()
{
  openLogs("scan");
  style();
  loadBand();
  ROOT::EnableImplicitMT(6);
  ROOT::RDataFrame d0("DecayTree", k10);
  auto d = derive(ROOT::RDF::RNode(d0));
  auto sels = scanSelections();
  std::vector<ROOT::RDF::RResultPtr<TH1D>> hs;
  for (size_t i = 0; i < sels.size(); ++i) hs.push_back(d.Filter(sels[i].cut).Histo1D({Form("hs%zu", i), "", kNb, kLo, kHi}, "D0_mass"));
  FILE *tab = fopen("tables/scan_results.tex", "w");
  int nPos2 = 0, nFreeLim = 0;
  double zmax = -99;
  std::string zmaxSel;
  for (size_t i = 0; i < sels.size(); ++i)
  {
    TH1D *h = hs[i].GetPtr();
    Fix fx = fitFixed(h, mD0, 0.015, true);
    Free fr = fitFree((TH1D *) h->Clone(Form("hc%zu", i)));
    auto [sb, esb] = sideband(h, mD0, 0.015);
    double z = fx.eS > 0 ? fx.S / fx.eS : 0;
    if (z > 2) ++nPos2;
    if (fr.atLimit) ++nFreeLim;
    if (z > zmax) { zmax = z; zmaxSel = sels[i].key; }
    L("%-4s %-45s N=%9.0f | fixed: S=%8.1f +- %6.1f (z=%5.2f) chi2/ndf=%.2f | sideband: S=%8.1f +- %6.1f | free: S=%8.1f mean=%.4f sigma=%.4f %s\n",
      sels[i].key.c_str(), sels[i].label.c_str(), h->GetEntries(), fx.S, fx.eS, z, fx.chi2 / fx.ndf, sb, esb, fr.S, fr.mean, fr.sig,
      fr.atLimit ? "AT LIMIT" : "");
    fprintf(tab, "%s & %s & %s & $%.0f\\pm%.0f$ & %.1f & %.1f & $%.0f\\pm%.0f$ & %.3f / %.0f%s \\\\\n", sels[i].key.c_str(), sels[i].label.c_str(),
            big(h->GetEntries()).c_str(), fx.S, fx.eS, z, fx.chi2 / fx.ndf, sb, esb, fr.mean, 1000 * fr.sig, fr.atLimit ? "$^\\dagger$" : "");
  }
  fclose(tab);
  V("scanN", Form("%zu", sels.size()));
  V("scanPosTwo", Form("%d", nPos2));
  V("scanFreeLim", Form("%d", nFreeLim));
  V("scanZmax", f1(zmax));
  V("scanZmaxSel", zmaxSel);
  // figures: 8 selections per figure
  for (int pg = 0; pg < 2; ++pg)
  {
    TCanvas c("c", "c", 1600, 760);
    c.Divide(4, 2);
    for (int k = 0; k < 8; ++k)
    {
      int i = pg * 8 + k;
      c.cd(k + 1);
      TH1D *h = hs[i].GetPtr();
      Fix fx = fitFixed(h, mD0, 0.015, false);
      drawMass(h, sels[i].key, fx.S, fx.eS);
    }
    save(c, Form("scan_mass_%d", pg + 1));
  }
  closeLogs();
}

// ---------------- grid scan
void stepGrid()
{
  openLogs("grid");
  style();
  loadBand();
  ROOT::EnableImplicitMT(6);
  ROOT::RDataFrame d0("DecayTree", k10);
  auto d = derive(ROOT::RDF::RNode(d0));
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
            if (pid) c += " && kPID && piPID";
            if (cts) c += " && abs(cts)<0.8";
            cfgs.push_back({Form("%.3f & %.0f & %.0f & %s & %s", dira, ls, pt, pid ? "yes" : "--", cts ? "yes" : "--"), c});
          }
  std::vector<ROOT::RDF::RResultPtr<TH1D>> hs;
  for (size_t i = 0; i < cfgs.size(); ++i) hs.push_back(d.Filter(cfgs[i].cut).Histo1D({Form("g%zu", i), "", kNb, kLo, kHi}, "D0_mass"));
  V("gridN", Form("%zu", cfgs.size()));
  TCanvas c("c", "c", 1500, 520);
  c.Divide(3, 1);
  std::vector<TH1D *> hz;
  const double sigs[3] = {0.010, 0.015, 0.022};
  const char *sname[3] = {"Ten", "Fifteen", "TwentyTwo"};
  for (int is = 0; is < 3; ++is)
  {
    double sig = sigs[is];
    struct R
    {
      size_t i;
      double S, eS, N, B;
    };
    std::vector<R> res;
    for (size_t i = 0; i < cfgs.size(); ++i)
    {
      TH1D *h = hs[i].GetPtr();
      if (h->GetEntries() < 200) continue;
      Fix fx = fitFixed(h, mD0, sig);
      res.push_back({i, fx.S, fx.eS, h->GetEntries(), fx.B});
    }
    std::sort(res.begin(), res.end(), [](const R &a, const R &b) { return a.S / a.eS > b.S / b.eS; });
    auto *h = new TH1D(Form("hz%d", is), ";S / #sigma_{S};selections", 24, -6, 6);
    int np2 = 0, nn2 = 0;
    double s1 = 0, s2 = 0;
    for (auto &r : res)
    {
      double z = r.S / r.eS;
      h->Fill(z);
      np2 += z > 2;
      nn2 += z < -2;
      s1 += z;
      s2 += z * z;
      L("sigma=%2.0f MeV  %-60s N=%8.0f S=%8.1f +- %6.1f z=%5.2f\n", sig * 1000, cfgs[r.i].cut.c_str(), r.N, r.S, r.eS, z);
    }
    double n = res.size(), mean = s1 / n, rms = std::sqrt(s2 / n - mean * mean);
    L("sigma=%2.0f MeV: %0.f fitted selections, z>2: %d, z<-2: %d, mean z %.2f, rms z %.2f\n", sig * 1000, n, np2, nn2, mean, rms);
    V(std::string("gridFit") + sname[is], f0(n));
    V(std::string("gridPos") + sname[is], Form("%d", np2));
    V(std::string("gridNeg") + sname[is], Form("%d", nn2));
    V(std::string("gridMean") + sname[is], f2(mean));
    V(std::string("gridRms") + sname[is], f2(rms));
    V(std::string("gridExp") + sname[is], f1(n * 0.02275));
    if (is == 1)
    {
      FILE *t = fopen("tables/grid_top.tex", "w");
      for (int k = 0; k < 8; ++k)
        fprintf(t, "%s & %s & $%.0f\\pm%.0f$ & %.2f \\\\\n", cfgs[res[k].i].name.c_str(), big(res[k].N).c_str(), res[k].S, res[k].eS, res[k].S / res[k].eS);
      fprintf(t, "\\midrule\n");
      for (size_t k = res.size() - 5; k < res.size(); ++k)
        fprintf(t, "%s & %s & $%.0f\\pm%.0f$ & %.2f \\\\\n", cfgs[res[k].i].name.c_str(), big(res[k].N).c_str(), res[k].S, res[k].eS, res[k].S / res[k].eS);
      fclose(t);
      V("gridTopZ", f2(res[0].S / res[0].eS));
      V("gridTopS", f0(res[0].S));
      V("gridTopES", f0(res[0].eS));
      V("gridBotZ", f2(res.back().S / res.back().eS));
    }
    c.cd(is + 1);
    h->SetLineWidth(2);
    h->SetMaximum(h->GetMaximum() * 1.7);
    h->Draw("hist");
    TF1 g("g", "[0]*exp(-0.5*x*x)/sqrt(2*pi)", -6, 6);
    g.SetParameter(0, n * h->GetBinWidth(1));
    g.SetLineColor(kRed + 1);
    g.SetLineStyle(2);
    g.DrawCopy("same");
    TLatex t;
    t.SetNDC();
    t.SetTextSize(0.05);
    t.DrawLatex(0.18, 0.87, Form("#sigma_{D^{0}} = %.0f MeV: rms %.2f", sig * 1000, rms));
    t.DrawLatex(0.18, 0.81, Form("z>2: %d, z<-2: %d (exp. %.1f each)", np2, nn2, n * 0.02275));
  }
  save(c, "grid_pulls");
  closeLogs();
}

// ---------------- local significance vs trial mass, background model pol2/3/4
void stepTrial()
{
  openLogs("trial");
  style();
  ROOT::EnableImplicitMT(6);
  ROOT::RDataFrame d("DecayTree", k10);
  std::string cut = "D0_DIRA>0.995 && D0_decayLength/D0_decayLengthErr>2";
  auto hp = d.Filter(cut).Histo1D({"h", "", kNb, kLo, kHi}, "D0_mass");
  TH1D *h = hp.GetPtr();
  L("selection: %s, entries %.0f\n", cut.c_str(), h->GetEntries());
  V("trialN", big(h->GetEntries()));
  TGraph g[3];
  const char *on[3] = {"Two", "Three", "Four"};
  double sig = 0.015;
  for (int ord = 2; ord <= 4; ++ord)
  {
    int ip = 0;
    double zmin = 99, zmaxv = -99, mzmax = 0;
    for (double m = 1.74; m <= 1.9901; m += 0.005)
    {
      TF1 f("f", Form("[0]*exp(-0.5*((x-[1])/[2])^2)+pol%d(3)", ord), 1.66, 2.08);
      f.SetParameter(0, 0);
      f.FixParameter(1, m);
      f.FixParameter(2, sig);
      f.SetParameter(3, h->GetBinContent(50));
      for (int k = 4; k < 4 + ord; ++k) f.SetParameter(k, 0);
      h->Fit(&f, "QN0RL");
      h->Fit(&f, "QN0RL");
      double z = f.GetParError(0) > 0 ? f.GetParameter(0) / f.GetParError(0) : 0;
      g[ord - 2].SetPoint(ip++, m, z);
      zmin = std::min(zmin, z);
      if (z > zmaxv) { zmaxv = z; mzmax = m; }
      if (std::fabs(m - 1.865) < 0.0026)
      {
        L("pol%d: m=%.3f z=%.2f chi2/ndf=%.2f\n", ord, m, z, f.GetChisquare() / f.GetNDF());
        V(std::string("trialZ") + on[ord - 2], f1(z));
        V(std::string("trialChi") + on[ord - 2], f1(f.GetChisquare() / f.GetNDF()));
      }
    }
    L("pol%d: z range %.1f .. %.1f, maximum at %.3f\n", ord, zmin, zmaxv, mzmax);
    V(std::string("trialZmin") + on[ord - 2], f1(zmin));
    V(std::string("trialZmax") + on[ord - 2], f1(zmaxv));
    V(std::string("trialMmax") + on[ord - 2], f3(mzmax));
  }
  TCanvas c("c", "c", 1400, 560);
  c.Divide(2, 1);
  c.cd(1);
  drawMass(h, "DIRA>0.995, L/#sigma_{L}>2", 0, 0, false);
  c.cd(2);
  gPad->SetLeftMargin(0.17);
  int col[3] = {kBlack, kBlue + 1, kRed + 1};
  TLegend lg(0.6, 0.16, 0.94, 0.36);
  lg.SetBorderSize(0);
  lg.SetTextSize(0.045);
  for (int k = 0; k < 3; ++k)
  {
    g[k].SetLineColor(col[k]);
    g[k].SetLineWidth(2);
    g[k].SetLineStyle(k + 1);
    lg.AddEntry(&g[k], Form("pol%d background", k + 2), "l");
  }
  g[0].GetXaxis()->SetTitle("trial mass (GeV/c^{2})");
  g[0].GetYaxis()->SetTitle("S / #sigma_{S}  (gaussian #sigma = 15 MeV)");
  g[0].GetYaxis()->SetTitleOffset(1.3);
  g[0].GetYaxis()->SetRangeUser(-16, 14);
  g[0].Draw("AL");
  g[1].Draw("L");
  g[2].Draw("L");
  TLine l(mD0, -16, mD0, 14);
  l.SetLineStyle(3);
  l.Draw();
  TLine l0(1.74, 0, 1.99, 0);
  l0.SetLineColor(kGray);
  l0.Draw();
  lg.Draw();
  save(c, "trial_mass");
  closeLogs();
}

// ---------------- ratio pointing / non-pointing
void stepRatio()
{
  openLogs("ratio");
  style();
  ROOT::EnableImplicitMT(6);
  std::string ctrl = "D0_DIRA>0.95 && D0_DIRA<0.98";
  struct S
  {
    std::string key, label, sig, ctl;
  };
  auto sets = [&](bool full) {
    std::string p = full ? "1" : "D0_pT>1.5";
    std::vector<S> v = {
        {"a", "DIRA>0.99, L/#sigma>2", p + " && D0_DIRA>0.99 && ls>2", p + " && " + ctrl},
        {"b", "DIRA>0.995, L/#sigma>2", p + " && D0_DIRA>0.995 && ls>2", p + " && " + ctrl},
        {"c", "DIRA>0.998, L/#sigma>2", p + " && D0_DIRA>0.998 && ls>2", p + " && " + ctrl},
        {"d", "DIRA>0.995, L/#sigma>4", p + " && D0_DIRA>0.995 && ls>4", p + " && " + ctrl},
        {"e", "DIRA>0.995, L/#sigma>2, p_{T} 1.5-2.5", p + " && D0_DIRA>0.995 && ls>2 && D0_pT<2.5", p + " && D0_pT<2.5 && " + ctrl},
        {"f", "DIRA>0.995, L/#sigma>2, p_{T}>2.5", p + " && D0_DIRA>0.995 && ls>2 && D0_pT>2.5", p + " && D0_pT>2.5 && " + ctrl},
        {"g", "control: DIRA 0.98-0.99", p + " && D0_DIRA>0.98 && D0_DIRA<0.99", p + " && " + ctrl},
    };
    if (!full) v.push_back({"h", "p_{T}<1.5: DIRA>0.995, L/#sigma>2", "D0_pT<1.5 && D0_DIRA>0.995 && ls>2", "D0_pT<1.5 && " + ctrl});
    return v;
  };
  FILE *tab = fopen("tables/ratio_results.tex", "w");
  std::vector<std::string> zrows;
  std::vector<std::pair<std::string, TGraph *>> trials;
  for (auto [tag, files] : std::vector<std::pair<std::string, const char *>>{{"Ten", k10}, {"Full", kFull}})
  {
    bool full = tag == "Full";
    ROOT::RDataFrame d0("DecayTree", files);
    auto d = d0.Define("ls", "D0_decayLength/D0_decayLengthErr");
    auto v = sets(full);
    std::vector<ROOT::RDF::RResultPtr<TH1D>> hs, hc;
    for (auto &s : v)
    {
      hs.push_back(d.Filter(s.sig).Histo1D({Form("s%s%s", tag.c_str(), s.key.c_str()), "", 50, kLo, kHi}, "D0_mass"));
      hc.push_back(d.Filter(s.ctl).Histo1D({Form("c%s%s", tag.c_str(), s.key.c_str()), "", 50, kLo, kHi}, "D0_mass"));
    }
    TCanvas c("c", "c", 1500, full ? 1000 : 1000);
    c.Divide(4, 2);
    fprintf(tab, "\\midrule\n\\multicolumn{9}{l}{\\textit{%s}}\\\\\n", full ? "full dataset, skim with $p_T>1.5$" : "10\\% sample");
    for (size_t i = 0; i < v.size(); ++i)
    {
      auto *r = new TH1D(Form("r%s%zu", tag.c_str(), i), "", 50, kLo, kHi);
      for (int b = 1; b <= 50; ++b)
      {
        double ns = hs[i]->GetBinContent(b), nc = hc[i]->GetBinContent(b);
        if (ns <= 0 || nc <= 0) continue;
        r->SetBinContent(b, ns / nc);
        r->SetBinError(b, ns / nc * std::sqrt(1. / ns + 1. / nc));
      }
      r->Scale(r->GetNbinsX() / r->Integral());
      // (1) pol2 only, (2) + fixed gaussian at the D0 mass (sigma 15 MeV, amplitude unbounded), (3) + free broad gaussian
      TF1 p2("p2", "pol2", 1.65, 2.08);
      r->Fit(&p2, "QR0");
      TF1 fx("fx", "[0]*exp(-0.5*((x-[1])/[2])^2)+pol2(3)", 1.65, 2.08);
      fx.SetParameters(0, mD0, 0.015, p2.GetParameter(0), p2.GetParameter(1), p2.GetParameter(2));
      fx.FixParameter(1, mD0);
      fx.FixParameter(2, 0.015);
      r->Fit(&fx, "QR0");
      TF1 fb("fb", "[0]*exp(-0.5*((x-[1])/[2])^2)+pol2(3)", 1.65, 2.08);
      fb.SetParameters(0.05, 1.88, 0.05, p2.GetParameter(0), p2.GetParameter(1), p2.GetParameter(2));
      fb.SetParLimits(0, -1, 1);
      fb.SetParLimits(1, 1.75, 2.0);
      fb.SetParLimits(2, 0.005, 0.25);
      r->Fit(&fb, "QR0");
      TF1 f4("f4", "[0]*exp(-0.5*((x-[1])/[2])^2)+pol4(3)", 1.65, 2.08);
      f4.SetParameters(0, mD0, 0.015, p2.GetParameter(0), p2.GetParameter(1), p2.GetParameter(2), 0, 0);
      f4.FixParameter(1, mD0);
      f4.FixParameter(2, 0.015);
      r->Fit(&f4, "QR0");
      TF1 fnb("fnb", "[0]*exp(-0.5*((x-[1])/[2])^2)+[3]*exp(-0.5*((x-[4])/[5])^2)+pol2(6)", 1.65, 2.08);
      fnb.SetParameters(0, mD0, 0.015, fb.GetParameter(0), fb.GetParameter(1), fb.GetParameter(2), fb.GetParameter(3), fb.GetParameter(4), fb.GetParameter(5));
      fnb.FixParameter(1, mD0);
      fnb.FixParameter(2, 0.015);
      fnb.SetParLimits(4, 1.75, 2.0);
      fnb.SetParLimits(5, 0.04, 0.25);
      r->Fit(&fnb, "QR0");
      double z4 = f4.GetParError(0) > 0 ? f4.GetParameter(0) / f4.GetParError(0) : 0;
      double znb = fnb.GetParError(0) > 0 ? fnb.GetParameter(0) / fnb.GetParError(0) : 0;
      L("[%s] %s narrow D0 gaussian (sigma 15 MeV fixed): pol2 z=%.2f chi2/ndf=%.1f/%d | pol4 z=%.2f chi2/ndf=%.1f/%d | broad gaus+pol2 z=%.2f chi2/ndf=%.1f/%d\n",
        tag.c_str(), v[i].key.c_str(), fx.GetParError(0) > 0 ? fx.GetParameter(0) / fx.GetParError(0) : 0., fx.GetChisquare(), fx.GetNDF(), z4,
        f4.GetChisquare(), f4.GetNDF(), znb, fnb.GetChisquare(), fnb.GetNDF());
      std::string kz = std::string("ratio") + tag + (char) toupper(v[i].key[0]);
      V(kz + "Zfour", f1(z4));
      V(kz + "Znb", f1(znb));
      V(kz + "ChiTwo", Form("%.0f/%d", fx.GetChisquare(), fx.GetNDF()));
      V(kz + "ChiFour", Form("%.0f/%d", f4.GetChisquare(), f4.GetNDF()));
      zrows.push_back(Form("%s & %s & %.1f & %.0f/%d & %.1f & %.0f/%d & %.1f & %.0f/%d \\\\\n", v[i].key.c_str(), tag == "Full" ? "full" : "10\\%",
                           fx.GetParError(0) > 0 ? fx.GetParameter(0) / fx.GetParError(0) : 0., fx.GetChisquare(), fx.GetNDF(), z4, f4.GetChisquare(),
                           f4.GetNDF(), znb, fnb.GetChisquare(), fnb.GetNDF()));
      // trial-mass scan (narrow gaussian + pol4) for the full sample
      if (full && (v[i].key == "a" || v[i].key == "b" || v[i].key == "e" || v[i].key == "g"))
      {
        auto *g = new TGraph();
        for (double m = 1.74; m <= 1.9901; m += 0.005)
        {
          TF1 ft("ft", "[0]*exp(-0.5*((x-[1])/[2])^2)+pol4(3)", 1.65, 2.08);
          ft.SetParameters(0, m, 0.015, f4.GetParameter(3), f4.GetParameter(4), f4.GetParameter(5), f4.GetParameter(6), f4.GetParameter(7));
          ft.FixParameter(1, m);
          ft.FixParameter(2, 0.015);
          r->Fit(&ft, "QR0");
          g->SetPoint(g->GetN(), m, ft.GetParError(0) > 0 ? ft.GetParameter(0) / ft.GetParError(0) : 0);
        }
        double zmx = -99, mmx = 0;
        for (int k = 0; k < g->GetN(); ++k)
          if (g->GetY()[k] > zmx) { zmx = g->GetY()[k]; mmx = g->GetX()[k]; }
        L("[%s] %s trial-mass scan (pol4): max z %.2f at %.3f; at 1.865: %.2f\n", tag.c_str(), v[i].key.c_str(), zmx, mmx, g->Eval(1.865));
        V(kz + "TrialZmax", f1(zmx));
        V(kz + "TrialMmax", f3(mmx));
        trials.push_back({v[i].key + ": " + v[i].label, g});
      }
      bool lim = fb.GetParameter(2) > 0.2499 || fb.GetParameter(2) < 0.0051 || fb.GetParameter(1) < 1.7501 || fb.GetParameter(1) > 1.9999;
      double zfx = fx.GetParError(0) > 0 ? fx.GetParameter(0) / fx.GetParError(0) : 0;
      double dchi = p2.GetChisquare() - fb.GetChisquare();
      L("[%s] %s %-40s Nsig %9.0f Nctl %9.0f | pol2 chi2/ndf %.1f/%d | fixed D0 gaus: amp %.4f+-%.4f (z=%.2f) | broad: amp %.4f+-%.4f mean %.4f+-%.4f "
        "sigma %.4f+-%.4f chi2 %.1f/%d dchi2 %.1f %s\n",
        tag.c_str(), v[i].key.c_str(), v[i].label.c_str(), hs[i]->GetEntries(), hc[i]->GetEntries(), p2.GetChisquare(), p2.GetNDF(), fx.GetParameter(0),
        fx.GetParError(0), zfx, fb.GetParameter(0), fb.GetParError(0), fb.GetParameter(1), fb.GetParError(1), fb.GetParameter(2), fb.GetParError(2),
        fb.GetChisquare(), fb.GetNDF(), dchi, lim ? "AT LIMIT" : "");
      std::string lab = v[i].label;
      for (auto [a, b] : std::vector<std::pair<std::string, std::string>>{{"#sigma", "$\\sigma$"}, {"p_{T}", "$p_T$"}})
        for (size_t pos; (pos = lab.find(a)) != std::string::npos;) lab.replace(pos, a.size(), b);
      fprintf(tab, "%s & %s & %s & %.0f/%d & $%.3f\\pm%.3f$ & $%.3f\\pm%.3f$ & $%.3f\\pm%.3f$ & $%.0f\\pm%.0f$%s & %.0f \\\\\n", v[i].key.c_str(), lab.c_str(),
              big(hs[i]->GetEntries()).c_str(), p2.GetChisquare(), p2.GetNDF(), fx.GetParameter(0), fx.GetParError(0), fb.GetParameter(0), fb.GetParError(0),
              fb.GetParameter(1), fb.GetParError(1), 1000 * fb.GetParameter(2), 1000 * fb.GetParError(2), lim ? "$^\\dagger$" : "", dchi);
      std::string kk = std::string("ratio") + tag + (char) toupper(v[i].key[0]);
      V(kk + "Amp", f3(fb.GetParameter(0)));
      V(kk + "eAmp", f3(fb.GetParError(0)));
      V(kk + "Mean", f3(fb.GetParameter(1)));
      V(kk + "Sig", f0(1000 * fb.GetParameter(2)));
      V(kk + "eSig", f0(1000 * fb.GetParError(2)));
      V(kk + "Dchi", f0(dchi));
      V(kk + "Zfix", f1(zfx));
      c.cd(i + 1);
      r->SetMarkerStyle(20);
      r->SetMarkerSize(0.6);
      r->GetXaxis()->SetTitle("m(K#pi) (GeV/c^{2})");
      r->GetXaxis()->SetNdivisions(505);
      r->GetYaxis()->SetTitle("ratio (normalised)");
      r->GetYaxis()->SetTitleOffset(1.45);
      gPad->SetLeftMargin(0.19);
      r->SetMinimum(r->GetMinimum() * 0.9);
      r->SetMaximum(r->GetMaximum() * 1.12);
      r->Draw("e");
      p2.SetLineColor(kGray + 2);
      p2.SetLineStyle(2);
      p2.DrawCopy("same");
      fb.SetLineColor(kRed + 1);
      fb.DrawCopy("same");
      TLine l(mD0, r->GetMinimum(), mD0, r->GetMaximum());
      l.SetLineStyle(3);
      l.DrawClone();
      TLatex t;
      t.SetNDC();
      t.SetTextSize(0.05);
      t.DrawLatex(0.18, 0.87, Form("(%s) %s", v[i].key.c_str(), v[i].label.c_str()));
      t.DrawLatex(0.18, 0.80, Form("#mu=%.3f #sigma=%.0f MeV #Delta#chi^{2}=%.0f", fb.GetParameter(1), 1000 * fb.GetParameter(2), dchi));
    }
    save(c, full ? "ratio_full" : "ratio_ten");
  }
  fclose(tab);
  FILE *tz = fopen("tables/ratio_narrow.tex", "w");
  for (auto &r : zrows) fputs(r.c_str(), tz);
  fclose(tz);
  TCanvas c2("c2", "c2", 900, 600);
  int col[4] = {kBlack, kBlue + 1, kGreen + 2, kRed + 1};
  TLegend lg(0.16, 0.68, 0.62, 0.92);
  lg.SetBorderSize(0);
  lg.SetTextSize(0.04);
  for (size_t k = 0; k < trials.size(); ++k)
  {
    auto *g = trials[k].second;
    g->SetLineColor(col[k % 4]);
    g->SetLineWidth(2);
    g->SetLineStyle(k + 1);
    if (k == 0)
    {
      g->GetXaxis()->SetTitle("trial mass (GeV/c^{2})");
      g->GetYaxis()->SetTitle("S / #sigma_{S}  (#sigma = 15 MeV, pol4)");
      g->GetYaxis()->SetRangeUser(-6, 11);
      g->Draw("AL");
    }
    else
      g->Draw("L");
    std::string lab = trials[k].first;
    lg.AddEntry(g, lab.c_str(), "l");
  }
  TLine l(mD0, -6, mD0, 11);
  l.SetLineStyle(3);
  l.Draw();
  TLine l0(1.74, 0, 1.99, 0);
  l0.SetLineColor(kGray);
  l0.Draw();
  lg.Draw();
  save(c2, "ratio_trial_full");
  closeLogs();
}

void d0cutscan(const char *step)
{
  gROOT->SetBatch(true);
  std::string s = step;
  if (s == "inputs") stepInputs();
  else if (s == "scan") stepScan();
  else if (s == "grid") stepGrid();
  else if (s == "trial") stepTrial();
  else if (s == "ratio") stepRatio();
  else printf("unknown step %s\n", step);
}
