// hump.C : origin of the broad enhancement in the pointing / non-pointing K pi mass ratio (full pT>1.5 skim)
// usage (note root): root -l -b -q 'scripts/hump.C+'   -> logs/hump.txt, logs/hump_values.txt, figures/hump_*.{pdf,png}
// Tests: (1) pointing probability vs opening angle; (2) mass ratio in opening-angle slices; (3) ratio predicted by
// weighting the control with P(pointing | opening angle, D0 pT) -> double ratio; (4) injection of a narrow peak to show
// the weighting keeps a narrow signal; (5) the same ratios for the swapped (K<->pi) and pi pi mass hypotheses.
#include <Math/GenVector/VectorUtil.h>
#include <Math/Vector4D.h>
#include <ROOT/RDataFrame.hxx>
#include <TCanvas.h>
#include <TF1.h>
#include <TGraphErrors.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
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
  void save(TCanvas &c, const char *n)
  {
    c.SaveAs(Form("figures/%s.pdf", n));
    c.SaveAs(Form("figures/%s.png", n));
  }
  TH1D *hratio(const TH1D *a, const TH1D *b, const char *name)
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
    return r;
  }
  struct NarrowFit
  {
    double z, chi2;
    int ndf;
    double broadDchi;
  };
  // narrow gaussian (sigma 15 MeV, fixed at m) + pol2, and pol2 vs pol2 + free broad gaussian
  NarrowFit fitRatio(TH1D *r, double m = mD0, double lo = 1.65, double hi = 2.08)
  {
    TF1 p2("p2", "pol2", lo, hi);
    r->Fit(&p2, "QR0");
    TF1 f("f", "[0]*exp(-0.5*((x-[1])/[2])^2)+pol2(3)", lo, hi);
    f.SetParameters(0, m, 0.015, p2.GetParameter(0), p2.GetParameter(1), p2.GetParameter(2));
    f.FixParameter(1, m);
    f.FixParameter(2, 0.015);
    r->Fit(&f, "QR0");
    TF1 b("b", "[0]*exp(-0.5*((x-[1])/[2])^2)+pol2(3)", lo, hi);
    b.SetParameters(0.05, 0.5 * (lo + hi), 0.08, p2.GetParameter(0), p2.GetParameter(1), p2.GetParameter(2));
    b.SetParLimits(1, lo, hi);
    b.SetParLimits(2, 0.005, 0.25);
    r->Fit(&b, "QR0");
    return {f.GetParError(0) > 0 ? f.GetParameter(0) / f.GetParError(0) : 0, f.GetChisquare(), f.GetNDF(), p2.GetChisquare() - b.GetChisquare()};
  }
}  // namespace

void hump()
{
  gL = fopen("logs/hump.txt", "w");
  gV = fopen("logs/hump_values.txt", "w");
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
  const double mpi = 0.13957, mK = 0.493677;
  ROOT::RDataFrame d0("DecayTree", kFull);
  auto d = d0.Define("ls", "D0_decayLength/D0_decayLengthErr")
               .Define("theta",
                       [](float a, float b, float c, float e, float f, float g) {
                         return (float) ROOT::Math::VectorUtil::Angle(V4(a, b, c, 0), V4(e, f, g, 0));
                       },
                       {"track_1_pT", "track_1_pseudorapidity", "track_1_phi", "track_2_pT", "track_2_pseudorapidity", "track_2_phi"})
               .Define("mswap", [=](float a, float b, float c, float e, float f, float g) { return (float) (V4(a, b, c, mK) + V4(e, f, g, mpi)).M(); },
                       {"track_1_pT", "track_1_pseudorapidity", "track_1_phi", "track_2_pT", "track_2_pseudorapidity", "track_2_phi"})
               .Define("mpipi", [=](float a, float b, float c, float e, float f, float g) { return (float) (V4(a, b, c, mpi) + V4(e, f, g, mpi)).M(); },
                       {"track_1_pT", "track_1_pseudorapidity", "track_1_phi", "track_2_pT", "track_2_pseudorapidity", "track_2_phi"})
               .Define("mcheck", [=](float a, float b, float c, float e, float f, float g) { return (float) (V4(a, b, c, mpi) + V4(e, f, g, mK)).M(); },
                       {"track_1_pT", "track_1_pseudorapidity", "track_1_phi", "track_2_pT", "track_2_pseudorapidity", "track_2_phi"});
  const std::string SIG = "D0_DIRA>0.995 && ls>2", CTL = "D0_DIRA>0.95 && D0_DIRA<0.98";
  auto ds = d.Filter(SIG), dc = d.Filter(CTL);

  // ---- pass 1: maps and basic histograms
  const double ptb[] = {1.5, 2.0, 2.5, 3.0, 4.0, 6.0, 20.0};
  auto mapS = ds.Histo2D({"mapS", "", 60, 0, 1.8, 6, ptb}, "theta", "D0_pT");
  auto mapC = dc.Histo2D({"mapC", "", 60, 0, 1.8, 6, ptb}, "theta", "D0_pT");
  auto thS = ds.Histo1D({"thS", "", 60, 0, 1.8}, "theta");
  auto thC = dc.Histo1D({"thC", "", 60, 0, 1.8}, "theta");
  auto mS = ds.Histo1D({"mS", "", 50, 1.6, 2.1}, "D0_mass"), mC = dc.Histo1D({"mC", "", 50, 1.6, 2.1}, "D0_mass");
  auto sS = ds.Histo1D({"sS", "", 50, 1.6, 2.1}, "mswap"), sC = dc.Histo1D({"sC", "", 50, 1.6, 2.1}, "mswap");
  auto pS = ds.Histo1D({"pS", "", 50, 1.2, 1.8}, "mpipi"), pC = dc.Histo1D({"pC", "", 50, 1.2, 1.8}, "mpipi");
  auto thm = dc.Profile1D({"thm", "", 50, 1.6, 2.1}, "D0_mass", "theta");
  auto diraTh = d.Filter("ls>2").Histo2D({"diraTh", "", 45, 0, 1.8, 50, 0.95, 1.0}, "theta", "D0_DIRA");
  auto lTh = dc.Profile1D({"lTh", "", 45, 0, 1.8}, "theta", "D0_decayLength");
  auto lThS = ds.Profile1D({"lThS", "", 45, 0, 1.8}, "theta", "D0_decayLength");
  auto mdiff = d.Define("dm", "abs(mcheck-D0_mass)").Max("dm");
  auto nS = ds.Count(), nC = dc.Count();
  // injection sample: control candidates with 1.855 < m < 1.875
  auto inj = dc.Filter("D0_mass>1.855 && D0_mass<1.875");
  auto injMap = inj.Histo2D({"injMap", "", 60, 0, 1.8, 6, ptb}, "theta", "D0_pT");
  auto injM = inj.Histo1D({"injM", "", 50, 1.6, 2.1}, "D0_mass");

  L("pointing (%s): %llu, control (%s): %llu; |recomputed mass - stored mass| max %.4f GeV\n", SIG.c_str(), *nS, CTL.c_str(), *nC, *mdiff);
  V("humpNsig", Form("%llu", *nS));
  V("humpNctl", Form("%llu", *nC));
  V("humpMassCheck", Form("%.1f", 1000 * *mdiff));

  // pointing probability (unnormalised ratio) vs theta, and theta slices
  auto *Pth = hratio(thS.GetPtr(), thC.GetPtr(), "Pth");
  Pth->Scale(1.0 * *nC / *nS);
  L("pointing / control (normalised to the total ratio) vs opening angle:\n");
  for (int i = 1; i <= Pth->GetNbinsX(); i += 6) L("  theta %.2f-%.2f rad: %.3f +- %.3f\n", Pth->GetXaxis()->GetBinLowEdge(i), Pth->GetXaxis()->GetBinUpEdge(i), Pth->GetBinContent(i), Pth->GetBinError(i));
  for (double t : {0.15, 0.45, 0.9, 1.4}) V(Form("humpP%s", t == 0.15 ? "A" : t == 0.45 ? "B" : t == 0.9 ? "C" : "D"), Form("%.2f", Pth->GetBinContent(Pth->FindBin(t))));
  L("mean opening angle of control candidates vs mass:\n");
  for (int i = 1; i <= thm->GetNbinsX(); i += 5) L("  m %.2f: <theta> %.3f rad\n", thm->GetBinCenter(i), thm->GetBinContent(i));
  V("humpThLo", Form("%.2f", thm->GetBinContent(thm->FindBin(1.65))));
  V("humpThHi", Form("%.2f", thm->GetBinContent(thm->FindBin(2.05))));
  L("mean decay length vs opening angle (control / pointing):\n");
  for (int i = 2; i <= lTh->GetNbinsX(); i += 8) L("  theta %.2f: L %.4f / %.4f cm\n", lTh->GetBinCenter(i), lTh->GetBinContent(i), lThS->GetBinContent(i));

  // ---- theta slices: quantiles of the control
  double q[5], pr[5] = {0.2, 0.4, 0.6, 0.8, 1.0};
  thC->GetQuantiles(5, q, pr);
  std::vector<double> edges = {0, q[0], q[1], q[2], q[3], 1.8};
  std::vector<ROOT::RDF::RResultPtr<TH1D>> slS, slC;
  for (int k = 0; k < 5; ++k)
  {
    std::string c = Form("theta>=%f && theta<%f", edges[k], edges[k + 1]);
    slS.push_back(ds.Filter(c).Histo1D({Form("slS%d", k), "", 25, 1.6, 2.1}, "D0_mass"));
    slC.push_back(dc.Filter(c).Histo1D({Form("slC%d", k), "", 25, 1.6, 2.1}, "D0_mass"));
  }

  // ---- pass 2: weight the control by P(theta, pT) = mapS/mapC (and with injection)
  auto makeW = [](const TH2D *s, const TH2D *c) {
    auto *w = (TH2D *) s->Clone(Form("w_%s", s->GetName()));
    w->Divide(c);
    return w;
  };
  TH2D *W = makeW(mapS.GetPtr(), mapC.GetPtr());
  auto *mapSi = (TH2D *) mapS->Clone("mapSi");
  const double Ninj = 3000;
  mapSi->Add(injMap.GetPtr(), Ninj / injMap->GetEntries());
  TH2D *Wi = makeW(mapSi, mapC.GetPtr());
  auto look = [](const TH2D *w) {
    return [w](float th, float pt) {
      int b = w->FindFixBin(th, std::min(pt, 19.9f));
      return w->GetBinContent(b);
    };
  };
  auto dw = dc.Define("w", look(W), {"theta", "D0_pT"}).Define("wi", look(Wi), {"theta", "D0_pT"});
  auto pM = dw.Histo1D({"pM", "", 50, 1.6, 2.1}, "D0_mass", "w");
  auto pMi = dw.Histo1D({"pMi", "", 50, 1.6, 2.1}, "D0_mass", "wi");
  auto pSw = dw.Histo1D({"pSw", "", 50, 1.6, 2.1}, "mswap", "w");
  auto pPP = dw.Histo1D({"pPP", "", 50, 1.2, 1.8}, "mpipi", "w");
  // the weights' statistical uncertainty is neglected in the predicted histograms (errors from sum of w^2)

  // ratios
  auto norm = [](TH1D *h) { h->Scale(h->GetNbinsX() / h->Integral()); };
  TH1D *R = hratio(mS.GetPtr(), mC.GetPtr(), "R");
  TH1D *Rp = hratio(pM.GetPtr(), mC.GetPtr(), "Rp");
  double sc = R->Integral();
  R->Scale(R->GetNbinsX() / sc);
  Rp->Scale(Rp->GetNbinsX() / Rp->Integral());
  TH1D *D = hratio(mS.GetPtr(), pM.GetPtr(), "D");
  norm(D);
  auto *mSi = (TH1D *) mS->Clone("mSi");
  mSi->Add(injM.GetPtr(), Ninj / injM->GetEntries());
  TH1D *Ri = hratio(mSi, mC.GetPtr(), "Ri");
  norm(Ri);
  TH1D *Di = hratio(mSi, pMi.GetPtr(), "Di");
  norm(Di);
  TH1D *RS = hratio(sS.GetPtr(), sC.GetPtr(), "RS");
  norm(RS);
  TH1D *DS = hratio(sS.GetPtr(), pSw.GetPtr(), "DS");
  norm(DS);
  TH1D *RP = hratio(pS.GetPtr(), pC.GetPtr(), "RP");
  norm(RP);
  TH1D *DP = hratio(pS.GetPtr(), pPP.GetPtr(), "DP");
  norm(DP);

  auto report = [&](const char *key, const char *what, TH1D *h, double m, double lo, double hi) {
    NarrowFit f = fitRatio(h, m, lo, hi);
    L("%-46s narrow (15 MeV at %.3f): z = %5.2f chi2/ndf %.1f/%d | broad-gaussian dchi2 over pol2 %.1f\n", what, m, f.z, f.chi2, f.ndf, f.broadDchi);
    V(std::string("hump") + key + "Z", Form("%.1f", f.z));
    V(std::string("hump") + key + "Chi", Form("%.0f/%d", f.chi2, f.ndf));
    V(std::string("hump") + key + "Dchi", Form("%.0f", f.broadDchi));
    return f;
  };
  report("R", "ratio pointing/control", R, mD0, 1.65, 2.08);
  report("Rp", "predicted ratio (theta,pT weighting)", Rp, mD0, 1.65, 2.08);
  report("D", "double ratio data/predicted", D, mD0, 1.65, 2.08);
  report("Ri", "ratio, 3000 injected at 1.865", Ri, mD0, 1.65, 2.08);
  report("Di", "double ratio, 3000 injected", Di, mD0, 1.65, 2.08);
  report("RS", "ratio, swapped mass", RS, mD0, 1.65, 2.08);
  report("DS", "double ratio, swapped mass", DS, mD0, 1.65, 2.08);
  report("RP", "ratio, pi pi mass (fixed 1.5)", RP, 1.5, 1.25, 1.78);
  report("DP", "double ratio, pi pi mass", DP, 1.5, 1.25, 1.78);
  // injected excess recovered in the double ratio: counts in 1.85-1.88 above a pol2 fitted outside the window
  auto excess = [&](TH1D *num, TH1D *den) {
    auto *dr = hratio(num, den, "tmp");
    TF1 f("fx", [](double *x, double *p) {
      if (x[0] > 1.85 && x[0] < 1.88) TF1::RejectPoint();
      return p[0] + p[1] * x[0] + p[2] * x[0] * x[0];
    }, 1.65, 2.08, 3);
    dr->Fit(&f, "QR0");
    double ex = 0;
    for (int i = num->FindBin(1.8501); i <= num->FindBin(1.8799); ++i) ex += num->GetBinContent(i) - den->GetBinContent(i) * f.Eval(num->GetBinCenter(i));
    delete dr;
    return ex;
  };
  double ex0 = excess(mS.GetPtr(), pM.GetPtr()), ex1 = excess(mSi, pMi.GetPtr());
  L("injection test: excess in 1.85-1.88 in the double ratio: %.0f without, %.0f with %.0f injected -> recovered %.0f (%.0f%%)\n", ex0, ex1, Ninj, ex1 - ex0,
    100 * (ex1 - ex0) / Ninj);
  V("humpInj", Form("%.0f", Ninj));
  V("humpInjRec", Form("%.0f", 100 * (ex1 - ex0) / Ninj));
  V("humpExcessData", Form("%.0f", ex0));
  // theta slices
  L("ratio in opening-angle slices (quintiles of the control):\n");
  std::vector<TH1D *> slR;
  for (int k = 0; k < 5; ++k)
  {
    TH1D *r = hratio(slS[k].GetPtr(), slC[k].GetPtr(), Form("slR%d", k));
    norm(r);
    TF1 p1("p1", "pol1", 1.6, 2.1), p2("p2", "pol2", 1.6, 2.1);
    r->Fit(&p1, "QR0");
    r->Fit(&p2, "QR0");
    NarrowFit f = fitRatio(r, mD0, 1.6, 2.1);
    L("  theta %.3f-%.3f: pol1 chi2/ndf %.1f/%d, pol2 %.1f/%d, narrow z %.2f, broad dchi2 %.1f, mass range of candidates %.2f-%.2f\n", edges[k], edges[k + 1],
      p1.GetChisquare(), p1.GetNDF(), p2.GetChisquare(), p2.GetNDF(), f.z, f.broadDchi, slC[k]->GetMean() - slC[k]->GetRMS(), slC[k]->GetMean() + slC[k]->GetRMS());
    V(Form("humpSl%cDchi", 'A' + k), Form("%.0f", f.broadDchi));
    V(Form("humpSl%cZ", 'A' + k), Form("%.1f", f.z));
    V(Form("humpSl%cLo", 'A' + k), Form("%.2f", edges[k]));
    V(Form("humpSl%cHi", 'A' + k), Form("%.2f", edges[k + 1]));
    slR.push_back(r);
  }

  // ---------------- figures
  {
    TCanvas c("c", "c", 1500, 520);
    c.Divide(3, 1);
    c.cd(1);
    gPad->SetLogz();
    gPad->SetRightMargin(0.12);
    diraTh->GetXaxis()->SetTitle("opening angle #theta_{12} (rad)");
    diraTh->GetYaxis()->SetTitle("DIRA");
    diraTh->DrawCopy("colz");
    c.cd(2);
    Pth->GetXaxis()->SetTitle("opening angle #theta_{12} (rad)");
    Pth->GetYaxis()->SetTitle("pointing / non-pointing (norm.)");
    Pth->SetMarkerStyle(20);
    Pth->SetMarkerSize(0.6);
    Pth->SetMinimum(0);
    Pth->DrawCopy("e");
    c.cd(3);
    thm->GetXaxis()->SetTitle("m(K#pi) (GeV/c^{2})");
    thm->GetYaxis()->SetTitle("#LT#theta_{12}#GT (rad), control");
    thm->SetMarkerStyle(20);
    thm->SetMarkerSize(0.6);
    thm->DrawCopy("e");
    save(c, "hump_theta");
  }
  {
    TCanvas c("c", "c", 1500, 520);
    c.Divide(3, 1);
    c.cd(1);
    R->SetMarkerStyle(20);
    R->SetMarkerSize(0.6);
    R->GetXaxis()->SetTitle("m(K#pi) (GeV/c^{2})");
    R->GetYaxis()->SetTitle("ratio (normalised)");
    R->DrawCopy("e");
    Rp->SetLineColor(kRed + 1);
    Rp->SetLineWidth(2);
    Rp->DrawCopy("hist same");
    TLegend lg(0.2, 0.75, 0.9, 0.92);
    lg.SetBorderSize(0);
    lg.SetTextSize(0.045);
    lg.AddEntry(R, "data: pointing / non-pointing", "pe");
    lg.AddEntry(Rp, "predicted from P(#theta_{12}, p_{T})", "l");
    lg.DrawClone();
    c.cd(2);
    for (TH1D *h : {D, Di})
    {
      h->SetMarkerStyle(h == D ? 20 : 24);
      h->SetMarkerSize(0.7);
      h->SetMarkerColor(h == D ? kBlack : kBlue + 1);
      h->SetLineColor(h == D ? kBlack : kBlue + 1);
    }
    Di->GetXaxis()->SetTitle("m(K#pi) (GeV/c^{2})");
    Di->GetYaxis()->SetTitle("data / predicted");
    Di->SetMinimum(0.96);
    Di->SetMaximum(1.06);
    Di->DrawCopy("e");
    D->DrawCopy("e same");
    TLegend lg2(0.2, 0.76, 0.9, 0.92);
    lg2.SetBorderSize(0);
    lg2.SetTextSize(0.045);
    lg2.AddEntry(D, "data", "pe");
    lg2.AddEntry(Di, Form("with %.0f injected at 1.855-1.875", Ninj), "pe");
    lg2.DrawClone();
    c.cd(3);
    for (TH1D *h : {RS, DS})
    {
      h->SetMarkerStyle(h == RS ? 20 : 24);
      h->SetMarkerSize(0.7);
      h->SetMarkerColor(h == RS ? kBlack : kRed + 1);
      h->SetLineColor(h == RS ? kBlack : kRed + 1);
    }
    RS->GetXaxis()->SetTitle("swapped m(K#pi) (GeV/c^{2})");
    RS->GetYaxis()->SetTitle("ratio (normalised)");
    RS->DrawCopy("e");
    DS->DrawCopy("e same");
    TLegend lg3(0.2, 0.76, 0.9, 0.92);
    lg3.SetBorderSize(0);
    lg3.SetTextSize(0.045);
    lg3.AddEntry(RS, "pointing / non-pointing", "pe");
    lg3.AddEntry(DS, "data / predicted", "pe");
    lg3.DrawClone();
    for (int k = 1; k <= 3; ++k)
    {
      c.cd(k);
      TLine l(mD0, gPad->GetUymin(), mD0, gPad->GetUymax());
      l.SetLineStyle(3);
      l.DrawClone();
    }
    save(c, "hump_predict");
  }
  {
    TCanvas c("c", "c", 1500, 520);
    c.Divide(3, 1);
    c.cd(1);
    TLegend lg(0.2, 0.62, 0.7, 0.93);
    lg.SetBorderSize(0);
    lg.SetTextSize(0.04);
    int col[5] = {kBlack, kBlue + 1, kGreen + 2, kOrange + 7, kRed + 1};
    for (int k = 0; k < 5; ++k)
    {
      slR[k]->SetLineColor(col[k]);
      slR[k]->SetMarkerColor(col[k]);
      slR[k]->SetMarkerStyle(20 + k);
      slR[k]->SetMarkerSize(0.7);
      slR[k]->GetXaxis()->SetTitle("m(K#pi) (GeV/c^{2})");
      slR[k]->GetYaxis()->SetTitle("ratio (normalised per slice)");
      slR[k]->SetMinimum(0.3);
      slR[k]->SetMaximum(2.2);
      slR[k]->DrawCopy(k ? "e same" : "e");
      lg.AddEntry(slR[k], Form("#theta_{12} %.2f-%.2f", edges[k], edges[k + 1]), "pe");
    }
    lg.DrawClone();
    c.cd(2);
    for (TH1D *h : {RP, DP})
    {
      h->SetMarkerStyle(h == RP ? 20 : 24);
      h->SetMarkerSize(0.7);
      h->SetMarkerColor(h == RP ? kBlack : kRed + 1);
      h->SetLineColor(h == RP ? kBlack : kRed + 1);
    }
    RP->GetXaxis()->SetTitle("m(#pi#pi) (GeV/c^{2})");
    RP->GetYaxis()->SetTitle("ratio (normalised)");
    RP->DrawCopy("e");
    DP->DrawCopy("e same");
    c.cd(3);
    lTh->SetMarkerStyle(24);
    lTh->SetMarkerSize(0.7);
    lThS->SetMarkerStyle(20);
    lThS->SetMarkerSize(0.7);
    lThS->SetMarkerColor(kRed + 1);
    lThS->SetLineColor(kRed + 1);
    lThS->GetXaxis()->SetTitle("opening angle #theta_{12} (rad)");
    lThS->GetYaxis()->SetTitle("#LT decay length #GT (cm)");
    lThS->SetMinimum(0);
    lThS->DrawCopy("e");
    lTh->DrawCopy("e same");
    TLegend lg4(0.45, 0.78, 0.95, 0.92);
    lg4.SetBorderSize(0);
    lg4.SetTextSize(0.045);
    lg4.AddEntry(lThS.GetPtr(), "pointing", "pe");
    lg4.AddEntry(lTh.GetPtr(), "non-pointing", "pe");
    lg4.DrawClone();
    save(c, "hump_slices");
  }
  fclose(gL);
  fclose(gV);
}
