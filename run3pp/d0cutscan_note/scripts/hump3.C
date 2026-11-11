// hump3.C : does the pointing probability of random (prompt) pairs depend on the daughter momentum sharing, and does
// weighting the non-pointing control in (opening angle, momentum asymmetry, pT) reproduce the hump?
// usage (note root): root -l -b -q 'scripts/hump3.C+'  -> logs/hump3.txt, logs/hump3_values.txt, figures/hump3_*.{pdf,png}
#include <Math/GenVector/VectorUtil.h>
#include <Math/Vector4D.h>
#include <ROOT/RDataFrame.hxx>
#include <TCanvas.h>
#include <TF1.h>
#include <TH1D.h>
#include <TH3D.h>
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
  TH1D *hratio(const TH1D *a, const TH1D *b, const char *name, bool nrm = true)
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
    if (nrm) r->Scale(r->GetNbinsX() / r->Integral());
    return r;
  }
  struct Fit
  {
    double chiP1, dchi, zn, ezn;
    int ndf;
  };
  Fit fitR(TH1D *r)
  {
    TF1 p1("p1", "pol1", 1.65, 2.08);
    r->Fit(&p1, "QR0");
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
    return {p1.GetChisquare(), p2.GetChisquare() - b.GetChisquare(), n.GetParameter(0), n.GetParError(0), p1.GetNDF()};
  }
}  // namespace

void hump3()
{
  gL = fopen("logs/hump3.txt", "w");
  gV = fopen("logs/hump3_values.txt", "w");
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
  std::vector<std::string> trk = {"track_1_pT", "track_1_pseudorapidity", "track_1_phi", "track_2_pT", "track_2_pseudorapidity", "track_2_phi"};
  ROOT::RDataFrame d0("DecayTree", kFull);
  auto d = d0.Define("ls", "D0_decayLength/D0_decayLengthErr")
               .Define("theta", [](float a, float b, float c, float x, float y, float z) { return (float) ROOT::Math::VectorUtil::Angle(V4(a, b, c, 0), V4(x, y, z, 0)); }, trk)
               .Define("asym", "float p1=track_1_pT*cosh(track_1_pseudorapidity), p2=track_2_pT*cosh(track_2_pseudorapidity); return (p1-p2)/(p1+p2);")
               .Filter("D0_decayLength<0.1", "prompt-like vertices (L < 1 mm)");
  const std::string SIG = "D0_DIRA>0.995 && ls>2", CTL = "D0_DIRA>0.95 && D0_DIRA<0.98";
  auto ds = d.Filter(SIG), dc = d.Filter(CTL);
  const double ptb[] = {1.5, 2.0, 2.5, 3.0, 4.0, 20.0};
  const int nb = 43;
  const double mlo = 1.6, mhi = 2.1 - 0.5 / 50 * 7;
  // 1D dependences
  auto aS = ds.Histo1D({"aS", "", 40, -1, 1}, "asym"), aC = dc.Histo1D({"aC", "", 40, -1, 1}, "asym");
  auto am = dc.Profile1D({"am", "", 43, mlo, mhi}, "D0_mass", "asym");
  // weight maps: (theta, asym, pT)
  struct Map
  {
    std::string key, label;
    int nth, nas;
  };
  std::vector<Map> maps = {{"T", "#theta_{12}, p_{T}", 30, 1}, {"A", "A, p_{T}", 1, 20}, {"TA", "#theta_{12}, A, p_{T}", 30, 20}, {"TAc", "#theta_{12}, A, p_{T} (coarse)", 10, 8}};
  std::vector<ROOT::RDF::RResultPtr<TH3D>> mS, mC, mI;
  auto inj = dc.Filter("D0_mass>1.855 && D0_mass<1.875");
  for (auto &m : maps)
  {
    TH3D model(("m" + m.key).c_str(), "", m.nth, 0, 1.8, m.nas, -1, 1, 5, 0, 1);
    model.GetZaxis()->Set(5, ptb);
    mS.push_back(ds.Histo3D(ROOT::RDF::TH3DModel(model), "theta", "asym", "D0_pT"));
    mC.push_back(dc.Histo3D(ROOT::RDF::TH3DModel(model), "theta", "asym", "D0_pT"));
    mI.push_back(inj.Histo3D(ROOT::RDF::TH3DModel(model), "theta", "asym", "D0_pT"));
  }
  auto hS = ds.Histo1D({"hS", "", nb, mlo, mhi}, "D0_mass"), hC = dc.Histo1D({"hC", "", nb, mlo, mhi}, "D0_mass");
  auto hI = inj.Histo1D({"hI", "", nb, mlo, mhi}, "D0_mass");
  auto nS = ds.Count(), nC = dc.Count();
  L("prompt-like (L<1 mm): pointing %llu, control %llu\n", *nS, *nC);
  V("hthreeNs", Form("%llu", *nS));
  TH1D *Ra = hratio(aS.GetPtr(), aC.GetPtr(), "Ra");
  L("pointing / non-pointing vs momentum asymmetry A=(p_pi-p_K)/(p_pi+p_K) (normalised):\n");
  for (int i = 2; i <= 40; i += 4) L("  A %+.2f: %.3f +- %.3f\n", Ra->GetBinCenter(i), Ra->GetBinContent(i), Ra->GetBinError(i));
  V("hthreeRaLo", Form("%.2f", Ra->GetBinContent(Ra->FindBin(-0.85))));
  V("hthreeRaMid", Form("%.2f", Ra->GetBinContent(Ra->FindBin(0.02))));
  V("hthreeRaHi", Form("%.2f", Ra->GetBinContent(Ra->FindBin(0.85))));
  L("mean asymmetry of non-pointing candidates vs mass:\n");
  for (int i = 1; i <= am->GetNbinsX(); i += 6) L("  m %.3f: <A> %+.3f\n", am->GetBinCenter(i), am->GetBinContent(i));
  V("hthreeAmLo", Form("%+.2f", am->GetBinContent(am->FindBin(1.65))));
  V("hthreeAmHi", Form("%+.2f", am->GetBinContent(am->FindBin(2.0))));

  const double Ninj = 3000;
  TH1D *R0 = hratio(hS.GetPtr(), hC.GetPtr(), "R0");
  Fit f0 = fitR(R0);
  L("data ratio: pol1 chi2/ndf %.0f/%d, broad dchi2 %.0f, narrow amp %.4f +- %.4f\n", f0.chiP1, f0.ndf, f0.dchi, f0.zn, f0.ezn);
  V("hthreeRChi", Form("%.0f/%d", f0.chiP1, f0.ndf));
  V("hthreeRDchi", Form("%.0f", f0.dchi));
  std::vector<TH1D *> dr, pr;
  for (size_t k = 0; k < maps.size(); ++k)
  {
    auto *W = (TH3D *) mS[k]->Clone(Form("W%zu", k));
    W->Divide(mC[k].GetPtr());
    auto *Si = (TH3D *) mS[k]->Clone(Form("Si%zu", k));
    Si->Add(mI[k].GetPtr(), Ninj / mI[k]->GetEntries());
    auto *Wi = (TH3D *) Si->Clone(Form("Wi%zu", k));
    Wi->Divide(mC[k].GetPtr());
    auto look = [](const TH3D *w) { return [w](float t, float a, float p) { return w->GetBinContent(w->FindFixBin(t, std::max(-0.9999f, std::min(a, 0.9999f)), std::min(p, 19.9f))); }; };
    auto dw = dc.Define("w", look(W), {"theta", "asym", "D0_pT"}).Define("wi", look(Wi), {"theta", "asym", "D0_pT"});
    auto P = dw.Histo1D({Form("P%zu", k), "", nb, mlo, mhi}, "D0_mass", "w");
    auto Pi = dw.Histo1D({Form("Pi%zu", k), "", nb, mlo, mhi}, "D0_mass", "wi");
    TH1D *D = hratio(hS.GetPtr(), P.GetPtr(), Form("D%zu", k));
    auto *hSi = (TH1D *) hS->Clone(Form("hSi%zu", k));
    hSi->Add(hI.GetPtr(), Ninj / hI->GetEntries());
    TH1D *Di = hratio(hSi, Pi.GetPtr(), Form("Di%zu", k));
    TH1D *Rp = hratio(P.GetPtr(), hC.GetPtr(), Form("Rp%zu", k));
    Fit fd = fitR(D), fi = fitR(Di);
    // excess counts in 1.845-1.885 (4 bins) over a pol2 fitted to the double ratio outside the window:
    // expected S_i = f(x_i) * P_i / k, with k the normalisation of the double ratio
    auto excess = [&](const TH1D *S, const TH1D *Pw, double &err) {
      TH1D *dd = hratio(S, Pw, "dd", false);
      double k = dd->GetNbinsX() / dd->Integral();
      dd->Scale(k);
      TF1 f("fx", [](double *x, double *p) {
        if (x[0] > 1.845 && x[0] < 1.885) TF1::RejectPoint();
        return p[0] + p[1] * x[0] + p[2] * x[0] * x[0];
      }, 1.65, 2.08, 3);
      dd->Fit(&f, "QR0");
      double ex = 0, n = 0;
      for (int i = S->GetXaxis()->FindFixBin(1.8451); i <= S->GetXaxis()->FindFixBin(1.8849); ++i)
      {
        ex += S->GetBinContent(i) - f.Eval(S->GetBinCenter(i)) * Pw->GetBinContent(i) / k;
        n += S->GetBinContent(i);
      }
      err = std::sqrt(n);
      delete dd;
      return ex;
    };
    double e0, e1;
    double ex0 = excess(hS.GetPtr(), P.GetPtr(), e0), ex1 = excess(hSi, Pi.GetPtr(), e1);
    double rec = ex1 - ex0;
    L("   excess in 1.845-1.885 (data): %.0f +- %.0f (stat., window counts)\n", ex0, e0);
    V("hthree" + maps[k].key + "Ex", Form("%.0f", ex0));
    V("hthree" + maps[k].key + "eEx", Form("%.0f", e0));
    L("weighting in (%s): predicted ratio pol1 chi2/ndf %.0f/%d | double ratio: pol1 chi2/ndf %.0f/%d, broad dchi2 %.0f, narrow z %.2f | injected %.0f -> recovered %.0f (%.0f%%), narrow z %.2f\n",
      maps[k].label.c_str(), fitR(Rp).chiP1, fd.ndf, fd.chiP1, fd.ndf, fd.dchi, fd.zn / fd.ezn, Ninj, rec, 100 * rec / Ninj, fi.zn / fi.ezn);
    V("hthree" + maps[k].key + "Chi", Form("%.0f/%d", fd.chiP1, fd.ndf));
    V("hthree" + maps[k].key + "Dchi", Form("%.0f", fd.dchi));
    V("hthree" + maps[k].key + "Z", Form("%.1f", fd.zn / fd.ezn));
    V("hthree" + maps[k].key + "Rec", Form("%.0f", 100 * rec / Ninj));
    V("hthree" + maps[k].key + "Zi", Form("%.1f", fi.zn / fi.ezn));
    dr.push_back(D);
    pr.push_back(Rp);
  }
  V("hthreeInj", Form("%.0f", Ninj));
  // figures
  TCanvas c("c", "c", 1500, 520);
  c.Divide(3, 1);
  c.cd(1);
  Ra->GetXaxis()->SetTitle("A = (p_{#pi} - p_{K}) / (p_{#pi} + p_{K})");
  Ra->GetYaxis()->SetTitle("pointing / non-pointing (norm.)");
  Ra->SetMarkerStyle(20);
  Ra->SetMarkerSize(0.7);
  Ra->SetMinimum(0);
  Ra->DrawCopy("e");
  c.cd(2);
  am->GetXaxis()->SetTitle("m(K#pi) (GeV/c^{2})");
  am->GetYaxis()->SetTitle("#LT A #GT, non-pointing");
  am->GetYaxis()->SetTitleOffset(1.5);
  gPad->SetLeftMargin(0.2);
  am->SetMarkerStyle(20);
  am->SetMarkerSize(0.7);
  am->DrawCopy("e");
  c.cd(3);
  R0->SetMarkerStyle(20);
  R0->SetMarkerSize(0.7);
  R0->GetXaxis()->SetTitle("m(K#pi) (GeV/c^{2})");
  R0->GetYaxis()->SetTitle("ratio (normalised)");
  R0->SetMinimum(0.8);
  R0->SetMaximum(1.25);
  R0->DrawCopy("e");
  int col[4] = {kBlue + 1, kGreen + 2, kRed + 1, kOrange + 7};
  TLegend lg(0.2, 0.72, 0.95, 0.93);
  lg.SetBorderSize(0);
  lg.SetTextSize(0.04);
  lg.AddEntry(R0, "data", "pe");
  for (size_t k = 0; k < pr.size(); ++k)
  {
    if (maps[k].key == "TAc") continue;
    pr[k]->SetLineColor(col[k]);
    pr[k]->SetLineWidth(2);
    pr[k]->DrawCopy("hist same");
    lg.AddEntry(pr[k], ("predicted, weights in " + maps[k].label).c_str(), "l");
  }
  lg.DrawClone();
  TLine l(mD0, 0.8, mD0, 1.25);
  l.SetLineStyle(3);
  l.DrawClone();
  c.SaveAs("figures/hump3_asym.pdf");
  c.SaveAs("figures/hump3_asym.png");
  TCanvas c2("c2", "c2", 900, 600);
  TLegend lg2(0.2, 0.75, 0.9, 0.93);
  lg2.SetBorderSize(0);
  lg2.SetTextSize(0.04);
  for (size_t k = 0; k < dr.size(); ++k)
  {
    dr[k]->SetMarkerStyle(20 + k);
    dr[k]->SetMarkerSize(0.8);
    dr[k]->SetMarkerColor(col[k]);
    dr[k]->SetLineColor(col[k]);
    dr[k]->GetXaxis()->SetTitle("m(K#pi) (GeV/c^{2})");
    dr[k]->GetYaxis()->SetTitle("data / predicted (normalised)");
    dr[k]->SetMinimum(0.9);
    dr[k]->SetMaximum(1.12);
    dr[k]->DrawCopy(k ? "e same" : "e");
    lg2.AddEntry(dr[k], ("weights in " + maps[k].label).c_str(), "pe");
  }
  lg2.DrawClone();
  TLine l2(mD0, 0.9, mD0, 1.12);
  l2.SetLineStyle(3);
  l2.DrawClone();
  c2.SaveAs("figures/hump3_double.pdf");
  c2.SaveAs("figures/hump3_double.png");
  fclose(gL);
  fclose(gV);
}
