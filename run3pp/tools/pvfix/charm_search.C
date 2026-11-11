// charm_search.C : signal search for 3-prong charm hadrons (D+ -> K- pi+ pi+, Ds+ -> phi(K+K-) pi+, Lambda_c+ -> p K- pi+,
// + c.c.) in the stock-PV charm production (output/charm10pct), with the method of the D0 analysis (pvfix_d0_v2.C):
//  - signal region: pointing (DIRA>0.995, L/sigma>2) or tight (DIRA>0.999, L/sigma>5); control: DIRA 0.95-0.98
//  - background prediction: control spectrum weighted with the pointing probability measured in bins of
//    (f_max = largest daughter momentum / sum of daughter momenta, mother pT)  [3-body analogue of the D0 asymmetry]
//  - window excess (+-20 MeV around the PDG mass) above a pol2 fitted to data/prediction outside the window(s);
//    narrow Gaussian (15 MeV fixed) z; injection test; local-z scan of the window centre over the fit range
//  - two mass versions: stored KFParticle mass ("st") and kappa-corrected mass recomputed from the daughters ("kc");
//    for kc the weights use the corrected pT and f_max as well
// Channel-specific: daughters' |DCA| as in the job (D+ 50 um, Ds 30 um, Lc 20 um); all daughters pT > 0.5;
// prompt-like L < 2 mm (D+, c tau = 312 um) or 1 mm (Ds, Lc); Ds: |m(KK) - 1.01946| < 10 MeV.
//   root -l -b -q 'charm_search.C+("<ch: Dplus|Ds|Lc>", "<file glob>", "<kappa csv>", "<out prefix>", "<extra cut>", "<label>")'
#include "kappa_apply.C"

#include <TCanvas.h>
#include <TF1.h>
#include <TH2D.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TStyle.h>

namespace
{
  FILE *cL, *cV;
  template <typename... A>
  void CL(const char *f, A... a)
  {
    printf(f, a...);
    fprintf(cL, f, a...);
  }
  void CV(const std::string &k, const std::string &v) { fprintf(cV, "%s %s\n", k.c_str(), v.c_str()); }
  struct Chan
  {
    std::string name, mom, t1, t2, t3, xlab;
    double mPdg, lo, hi, fitLo, fitHi, dcaMin, lMax;
    std::vector<double> excl;  // other window centres excluded from the background fit (e.g. D+ -> KKpi in the Ds channel)
  };
  Chan chanOf(const std::string &c)
  {
    if (c == "Dplus") return {"Dplus", "Dplus", "track_1", "track_2", "track_3", "m(K#pi#pi) (GeV/c^{2})", 1.86966, 1.70, 2.05, 1.72, 2.04, 0.005, 0.2, {}};
    if (c == "Ds") return {"Ds", "D_splus", "phi_track_1", "phi_track_2", "track_3", "m(KK#pi) (GeV/c^{2})", 1.96835, 1.75, 2.10, 1.77, 2.09, 0.003, 0.1, {1.86966}};
    return {"Lc", "Lambda_cplus", "track_1", "track_2", "track_3", "m(pK#pi) (GeV/c^{2})", 2.28646, 2.15, 2.42, 2.17, 2.41, 0.002, 0.1, {}};
  }
  Chan gC;
  const double kHalf = 0.020;
  TH1D *cratio(const TH1D *a, const TH1D *b, const char *name, bool nrm = true)
  {
    auto *r = (TH1D *) a->Clone(name);
    r->Reset();
    for (int i = 1; i <= a->GetNbinsX(); ++i)
    {
      double x = a->GetBinContent(i), y = b->GetBinContent(i), ex = a->GetBinError(i), ey = b->GetBinError(i);
      if (x <= 0 || y <= 0) continue;
      r->SetBinContent(i, x / y);
      r->SetBinError(i, x / y * std::sqrt(ex * ex / (x * x) + ey * ey / (y * y)));
    }
    if (nrm && r->Integral() > 0) r->Scale(r->GetNbinsX() / r->Integral());
    return r;
  }
  double gCentre;
  bool inExcl(double x)
  {
    if (std::fabs(x - gCentre) < kHalf) return true;
    for (double e : gC.excl)
      if (std::fabs(e - gCentre) > 1e-6 && std::fabs(x - e) < kHalf) return true;
    return false;
  }
  struct CRes
  {
    double ex, eex, z, dchi, chi2;
    int ndf;
  };
  CRes canalyse(const TH1D *S, const TH1D *P, const char *tag, double centre, bool full = true)
  {
    gCentre = centre;
    TH1D *dd = cratio(S, P, Form("dd_%s_%.3f", tag, centre), false);
    double k = dd->GetNbinsX() / dd->Integral();
    dd->Scale(k);
    TF1 f(Form("fx_%s", tag), [](double *x, double *p) {
      if (inExcl(x[0])) TF1::RejectPoint();
      return p[0] + p[1] * x[0] + p[2] * x[0] * x[0];
    }, gC.fitLo, gC.fitHi, 3);
    dd->Fit(&f, "QR0");
    double ex = 0, n = 0;
    for (int i = S->GetXaxis()->FindFixBin(centre - kHalf + 1e-4); i <= S->GetXaxis()->FindFixBin(centre + kHalf - 1e-4); ++i)
    {
      ex += S->GetBinContent(i) - f.Eval(S->GetBinCenter(i)) * P->GetBinContent(i) / k;
      n += S->GetBinContent(i);
    }
    CRes r{ex, std::sqrt(n), 0, 0, f.GetChisquare(), f.GetNDF()};
    if (!full) return r;
    TF1 p2("p2", "pol2", gC.fitLo, gC.fitHi);
    dd->Fit(&p2, "QR0");
    TF1 g("g", "[0]*exp(-0.5*((x-[1])/[2])^2)+pol2(3)", gC.fitLo, gC.fitHi);
    g.SetParameters(0, centre, 0.015, p2.GetParameter(0), p2.GetParameter(1), p2.GetParameter(2));
    g.FixParameter(1, centre);
    g.FixParameter(2, 0.015);
    dd->Fit(&g, "QR0");
    TF1 b("b", "[0]*exp(-0.5*((x-[1])/[2])^2)+pol2(3)", gC.fitLo, gC.fitHi);
    b.SetParameters(0.05, 0.5 * (gC.fitLo + gC.fitHi), 0.08, p2.GetParameter(0), p2.GetParameter(1), p2.GetParameter(2));
    b.SetParLimits(1, gC.fitLo + 0.05, gC.fitHi - 0.05);
    b.SetParLimits(2, 0.005, 0.25);
    dd->Fit(&b, "QR0");
    r.z = g.GetParError(0) > 0 ? g.GetParameter(0) / g.GetParError(0) : 0;
    r.dchi = p2.GetChisquare() - b.GetChisquare();
    return r;
  }
  float hypoM(int pdg)
  {
    int a = std::abs(pdg);
    return a == 321 ? 0.493677f : a == 2212 ? 0.938272f : 0.13957f;
  }
  float mass3(float p1, float e1, float f1, float m1, float p2, float e2, float f2, float m2, float p3, float e3, float f3, float m3)
  {
    ROOT::Math::PtEtaPhiMVector a(p1, e1, f1, m1), b(p2, e2, f2, m2), c(p3, e3, f3, m3);
    return (float) (a + b + c).M();
  }
  float pt3(float p1, float f1, float p2, float f2, float p3, float f3)
  {
    return (float) std::hypot(p1 * std::cos(f1) + p2 * std::cos(f2) + p3 * std::cos(f3), p1 * std::sin(f1) + p2 * std::sin(f2) + p3 * std::sin(f3));
  }
  float fmax3(float p1, float e1, float p2, float e2, float p3, float e3)
  {
    float a = p1 * std::cosh(e1), b = p2 * std::cosh(e2), c = p3 * std::cosh(e3);
    return std::max(a, std::max(b, c)) / (a + b + c);
  }
}  // namespace

void charm_search(const char *ch, const char *in, const char *csv, const char *prefix, const char *extra = "1", const char *extraLabel = "")
{
  gC = chanOf(ch);
  loadMap(csv);
  cL = fopen(Form("%s.txt", prefix), "w");
  cV = fopen(Form("%s_values.txt", prefix), "w");
  ROOT::EnableImplicitMT(8);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetTitleSize(0.05, "XYZ");
  gStyle->SetLabelSize(0.045, "XYZ");
  gStyle->SetPadLeftMargin(0.16);
  gStyle->SetPadBottomMargin(0.13);
  gStyle->SetPadRightMargin(0.04);
  gStyle->SetPadTopMargin(0.07);
  const Chan &C = gC;
  const std::string M = C.mom, t[3] = {C.t1, C.t2, C.t3};
  ROOT::RDataFrame d0("DecayTree", in);
  ROOT::RDF::RNode d = d0;
  for (int i = 0; i < 3; ++i)
  {
    std::string k = Form("k%d", i + 1), q = t[i] + "_PDG_ID";
    d = d.Define(k, [](int pdg, float p, float e, float f) { return kappaOf(pdg > 0 ? 1 : -1, p, e, f); }, {q, t[i] + "_pT", t[i] + "_pseudorapidity", t[i] + "_phi"})
            .Define(Form("pc%d", i + 1), k + "*" + t[i] + "_pT")
            .Define(Form("mh%d", i + 1), [](int p) { return hypoM(p); }, {q});
  }
  auto cols = [&](const char *p) {
    return std::vector<std::string>{std::string(p) + "1", t[0] + "_pseudorapidity", t[0] + "_phi", "mh1", std::string(p) + "2", t[1] + "_pseudorapidity", t[1] + "_phi", "mh2",
                                    std::string(p) + "3", t[2] + "_pseudorapidity", t[2] + "_phi", "mh3"};
  };
  d = d.Alias("ps1", t[0] + "_pT").Alias("ps2", t[1] + "_pT").Alias("ps3", t[2] + "_pT");
  d = d.Define("m_rc", mass3, cols("ps")).Define("m_kc", mass3, cols("pc")).Define("m_st", M + "_mass");
  d = d.Define("pt_st", M + "_pT")
          .Define("pt_kc", pt3, {"pc1", t[0] + "_phi", "pc2", t[1] + "_phi", "pc3", t[2] + "_phi"})
          .Define("fm_st", fmax3, {"ps1", t[0] + "_pseudorapidity", "ps2", t[1] + "_pseudorapidity", "ps3", t[2] + "_pseudorapidity"})
          .Define("fm_kc", fmax3, {"pc1", t[0] + "_pseudorapidity", "pc2", t[1] + "_pseudorapidity", "pc3", t[2] + "_pseudorapidity"});
  std::string base = Form("ps1>0.5 && ps2>0.5 && ps3>0.5 && abs(%s_PV_DCA)>%g && abs(%s_PV_DCA)>%g && abs(%s_PV_DCA)>%g && %s_decayLength<%g", t[0].c_str(), C.dcaMin, t[1].c_str(),
                          C.dcaMin, t[2].c_str(), C.dcaMin, M.c_str(), C.lMax);
  if (C.name == "Ds") base += " && abs(phi_mass-1.01946)<0.010";
  base = "(" + std::string(extra) + ") && " + base;
  auto nAll = d0.Count();
  auto db = d.Filter(base);
  auto nBase = db.Count();
  auto drc = db.Filter(M + "_decayLength<0.1").Define("x", "1e3f*(m_rc-m_st)").Stats("x");
  std::string ls = "(" + M + "_decayLength/" + M + "_decayLengthErr)";
  struct Ver
  {
    std::string key, m, pt, fm, label;
  };
  std::vector<Ver> vers = {{"st", "m_st", "pt_st", "fm_st", "stored momenta"}, {"kc", "m_kc", "pt_kc", "fm_kc", "#kappa-corrected"}};
  struct Rng
  {
    std::string key, label;
  };
  std::vector<Rng> rngs = {{"Hi", "p_{T} > 1.5"}, {"All", "all p_{T}"}};
  struct Sel
  {
    std::string key, cut, label;
  };
  std::vector<Sel> sels = {{"P", M + "_DIRA>0.995 && " + ls + ">2", "DIRA>0.995, L/#sigma>2"}, {"T", M + "_DIRA>0.999 && " + ls + ">5", "DIRA>0.999, L/#sigma>5"}};
  const double ptb[] = {0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 20.0};
  const int nb = (int) std::lround((C.hi - C.lo) / 0.010);
  struct Book
  {
    std::string tag, vk, rk, sk;
    ROOT::RDF::RResultPtr<TH2D> mS, mC, mI;
    ROOT::RDF::RResultPtr<TH1D> hS, hC, hI;
    ROOT::RDF::RNode dc;
    std::string m, pt, fm;
  };
  std::vector<Book> bk;
  for (auto &v : vers)
    for (auto &r : rngs)
      for (auto &s : sels)
      {
        std::string rc = r.key == "Hi" ? v.pt + ">1.5" : "1";
        auto ds = db.Filter(rc + " && " + s.cut);
        auto dc = db.Filter(rc + " && " + M + "_DIRA>0.95 && " + M + "_DIRA<0.98");
        auto inj = dc.Filter(Form("%s>%g && %s<%g", v.m.c_str(), C.mPdg - 0.010, v.m.c_str(), C.mPdg + 0.010));
        std::string tag = v.key + r.key + s.key;
        TH2D model(("m" + tag).c_str(), "", 20, 1. / 3, 1, 7, ptb);
        bk.push_back({tag, v.key, r.key, s.key, ds.Histo2D(ROOT::RDF::TH2DModel(model), v.fm, v.pt), dc.Histo2D(ROOT::RDF::TH2DModel(model), v.fm, v.pt),
                      inj.Histo2D(ROOT::RDF::TH2DModel(model), v.fm, v.pt), ds.Histo1D({("hS" + tag).c_str(), "", nb, C.lo, C.hi}, v.m),
                      dc.Histo1D({("hC" + tag).c_str(), "", nb, C.lo, C.hi}, v.m), inj.Histo1D({("hI" + tag).c_str(), "", nb, C.lo, C.hi}, v.m), dc, v.m, v.pt, v.fm});
      }
  CL("channel %s (%s), PDG mass %.5f, window +-%.0f MeV; extra selection: %s [%s]\n", C.name.c_str(), M.c_str(), C.mPdg, 1e3 * kHalf, extra, extraLabel);
  CL("candidates %llu, after the offline base selection %llu (%.2f%%)\n", *nAll, *nBase, 100. * *nBase / *nAll);
  CL("recomputed - stored mass (L < 1 mm): mean %.2f MeV, rms %.2f MeV\n", drc->GetMean(), drc->GetRMS());
  CV("N", Form("%llu", *nAll));
  CV("Nbase", Form("%llu", *nBase));
  std::vector<TH1D *> spec, dblR;
  for (auto &b : bk)
  {
    double nS = b.hS->GetEntries(), nC = b.hC->GetEntries();
    const double Ninj = std::max(50.0, std::round(0.002 * nS));
    auto *W = (TH2D *) b.mS->Clone(("W" + b.tag).c_str());
    W->Divide(b.mC.GetPtr());
    auto *Si = (TH2D *) b.mS->Clone(("Si" + b.tag).c_str());
    Si->Add(b.mI.GetPtr(), Ninj / std::max(1.0, b.mI->GetEntries()));
    auto *Wi = (TH2D *) Si->Clone(("Wi" + b.tag).c_str());
    Wi->Divide(b.mC.GetPtr());
    auto look = [](const TH2D *w) { return [w](float a, float p) { return w->GetBinContent(w->FindFixBin(std::max(0.3334f, std::min(a, 0.9999f)), std::min(p, 19.9f))); }; };
    auto P = b.dc.Define("w", look(W), {b.fm, b.pt}).Histo1D({("P" + b.tag).c_str(), "", nb, C.lo, C.hi}, b.m, "w");
    auto Pi = b.dc.Define("wi", look(Wi), {b.fm, b.pt}).Histo1D({("Pi" + b.tag).c_str(), "", nb, C.lo, C.hi}, b.m, "wi");
    auto *hSi = (TH1D *) b.hS->Clone(("hSi" + b.tag).c_str());
    hSi->Add(b.hI.GetPtr(), Ninj / std::max(1.0, b.hI->GetEntries()));
    CRes r0 = canalyse(b.hS.GetPtr(), P.GetPtr(), b.tag.c_str(), C.mPdg), r1 = canalyse(hSi, Pi.GetPtr(), (b.tag + "i").c_str(), C.mPdg);
    double rec = r1.ex - r0.ex;
    // local-z scan of the window centre
    std::string scan;
    double zmaxScan = 0, cmax = 0;
    for (double c = C.fitLo + 0.03; c <= C.fitHi - 0.03 + 1e-6; c += 0.010)
    {
      CRes rs = canalyse(b.hS.GetPtr(), P.GetPtr(), (b.tag + "s").c_str(), c, false);
      double z = rs.eex > 0 ? rs.ex / rs.eex : 0;
      scan += Form(" %.3f:%.1f", c, z);
      if (std::fabs(z) > std::fabs(zmaxScan)) zmaxScan = z, cmax = c;
    }
    CL("%-6s %-16s %-10s %-18s Nsig %8.0f Nctl %8.0f | pol2 chi2/ndf %4.0f/%d, broad dchi2 %4.0f | excess %6.0f +- %4.0f (%.1f sigma), narrow z %5.2f | injected %.0f -> recovered %.0f%%\n",
       b.tag.c_str(), b.vk == "st" ? "stored momenta" : "kappa-corrected", b.rk == "Hi" ? "pT>1.5" : "all pT", b.sk == "P" ? "DIRA>0.995,L/s>2" : "DIRA>0.999,L/s>5", nS, nC, r0.chi2,
       r0.ndf, r0.dchi, r0.ex, r0.eex, r0.ex / r0.eex, r0.z, Ninj, 100 * rec / Ninj);
    CL("SCAN %s largest |z| %.1f at %.3f;%s\n", b.tag.c_str(), zmaxScan, cmax, scan.c_str());
    for (double e : C.excl)
    {
      CRes re = canalyse(b.hS.GetPtr(), P.GetPtr(), (b.tag + "e").c_str(), e);
      CL("ALT %s window at %.5f (D+ -> KKpi): excess %6.0f +- %4.0f (%.1f sigma), narrow z %5.2f\n", b.tag.c_str(), e, re.ex, re.eex, re.ex / re.eex, re.z);
      CV("c" + b.tag + "AltZ", Form("%.1f", re.ex / re.eex));
    }
    CV("c" + b.tag + "Ns", Form("%.0f", nS));
    CV("c" + b.tag + "Nc", Form("%.0f", nC));
    CV("c" + b.tag + "Ex", Form("%.0f", r0.ex));
    CV("c" + b.tag + "eEx", Form("%.0f", r0.eex));
    CV("c" + b.tag + "Z", Form("%.1f", r0.ex / r0.eex));
    CV("c" + b.tag + "Nz", Form("%.2f", r0.z));
    CV("c" + b.tag + "Chi", Form("%.0f/%d", r0.chi2, r0.ndf));
    CV("c" + b.tag + "Dchi", Form("%.0f", r0.dchi));
    CV("c" + b.tag + "Rec", Form("%.0f", 100 * rec / Ninj));
    CV("c" + b.tag + "ScanMax", Form("%.1f@%.3f", zmaxScan, cmax));
    spec.push_back((TH1D *) b.hS->Clone(("sp" + b.tag).c_str()));
    dblR.push_back(cratio(b.hS.GetPtr(), P.GetPtr(), ("D" + b.tag).c_str()));
  }
  int col[2] = {kBlack, kRed + 1};
  for (auto &r : rngs)
    for (auto &s : sels)
    {
      TCanvas c("c", "c", 1000, 500);
      c.Divide(2, 1);
      TLegend lg(0.2, 0.75, 0.7, 0.92);
      lg.SetBorderSize(0);
      lg.SetTextSize(0.042);
      double ymax = 0;
      std::vector<size_t> idx;
      for (int iv = 0; iv < 2; ++iv)
        for (size_t i = 0; i < bk.size(); ++i)
          if (bk[i].vk == vers[iv].key && bk[i].rk == r.key && bk[i].sk == s.key) idx.push_back(i), ymax = std::max(ymax, spec[i]->GetMaximum());
      for (int iv = 0; iv < 2; ++iv)
      {
        size_t i = idx[iv];
        for (int p = 0; p < 2; ++p)
        {
          TH1D *x = p == 0 ? spec[i] : dblR[i];
          x->SetLineColor(col[iv]);
          x->SetMarkerColor(col[iv]);
          x->SetMarkerStyle(20 + iv);
          x->SetMarkerSize(0.6);
          x->GetXaxis()->SetTitle(C.xlab.c_str());
          x->GetYaxis()->SetTitle(p == 0 ? "candidates / 10 MeV" : "data / weighted prediction (norm.)");
          if (p == 0) x->SetMinimum(0), x->SetMaximum(1.3 * ymax);
          else x->SetMinimum(0.8), x->SetMaximum(1.2);
          c.cd(p + 1);
          x->DrawCopy(iv ? "e same" : "e");
        }
        lg.AddEntry(spec[i], vers[iv].label.c_str(), "pe");
      }
      for (int p = 1; p <= 2; ++p)
      {
        c.cd(p);
        gPad->Update();
        std::vector<double> ms = {C.mPdg};
        for (double e : C.excl) ms.push_back(e);
        for (double m : ms)
        {
          TLine l(m, gPad->GetUymin(), m, gPad->GetUymax());
          l.SetLineStyle(3);
          l.DrawClone();
        }
      }
      c.cd(1);
      lg.DrawClone();
      TLatex tx;
      tx.SetNDC();
      tx.SetTextSize(0.04);
      tx.DrawLatex(0.17, 0.95, Form("%s: %s, %s, L < %.0f mm %s", C.name.c_str(), r.label.c_str(), s.label.c_str(), 10 * C.lMax, extraLabel));
      c.SaveAs(Form("%s_%s%s.png", prefix, r.key.c_str(), s.key.c_str()));
    }
  fclose(cL);
  fclose(cV);
}
