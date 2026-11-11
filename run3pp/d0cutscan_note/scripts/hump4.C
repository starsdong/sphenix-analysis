// hump4.C : D0-window excess with the (A, pT)-weighted non-pointing background, in bins of PV / bunch-crossing multiplicity
// usage (note root): root -l -b -q 'scripts/hump4.C+'  -> logs/hump4.txt, logs/hump4_values.txt, figures/hump4_mult.{pdf,png}
#include <Math/Vector4D.h>
#include <ROOT/RDataFrame.hxx>
#include <TCanvas.h>
#include <TF1.h>
#include <TGraphErrors.h>
#include <TH1D.h>
#include <TH2D.h>
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
  const int nb = 43;
  const double mlo = 1.6, mhi = 2.1 - 0.5 / 50 * 7;
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
      double x = a->GetBinContent(i), y = b->GetBinContent(i), ex = a->GetBinError(i), ey = b->GetBinError(i);
      if (x <= 0 || y <= 0) continue;
      r->SetBinContent(i, x / y);
      r->SetBinError(i, x / y * std::sqrt(ex * ex / (x * x) + ey * ey / (y * y)));
    }
    return r;
  }
  // excess in 1.845-1.885 above pol2 fitted to data/predicted outside the window; z of a 15 MeV gaussian at m_D0
  struct Res
  {
    double ex, eex, z, dchi;
  };
  Res analyse(TH1D *S, TH1D *P)
  {
    TH1D *dd = hratio(S, P, Form("dd_%s", S->GetName()));
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
      ex += S->GetBinContent(i) - f.Eval(S->GetBinCenter(i)) * P->GetBinContent(i) / k;
      n += S->GetBinContent(i);
    }
    TF1 p2("p2", "pol2", 1.65, 2.08);
    dd->Fit(&p2, "QR0");
    TF1 g("g", "[0]*exp(-0.5*((x-[1])/[2])^2)+pol2(3)", 1.65, 2.08);
    g.SetParameters(0, mD0, 0.015, p2.GetParameter(0), p2.GetParameter(1), p2.GetParameter(2));
    g.FixParameter(1, mD0);
    g.FixParameter(2, 0.015);
    dd->Fit(&g, "QR0");
    TF1 b("b", "[0]*exp(-0.5*((x-[1])/[2])^2)+pol2(3)", 1.65, 2.08);
    b.SetParameters(0.05, 1.88, 0.08, p2.GetParameter(0), p2.GetParameter(1), p2.GetParameter(2));
    b.SetParLimits(1, 1.75, 2.0);
    b.SetParLimits(2, 0.005, 0.25);
    dd->Fit(&b, "QR0");
    return {ex, std::sqrt(n), g.GetParError(0) > 0 ? g.GetParameter(0) / g.GetParError(0) : 0, p2.GetChisquare() - b.GetChisquare()};
  }
}  // namespace

void hump4()
{
  gL = fopen("logs/hump4.txt", "w");
  gV = fopen("logs/hump4_values.txt", "w");
  ROOT::EnableImplicitMT(6);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetTitleSize(0.055, "XYZ");
  gStyle->SetLabelSize(0.05, "XYZ");
  gStyle->SetPadLeftMargin(0.17);
  gStyle->SetPadBottomMargin(0.14);
  gStyle->SetPadRightMargin(0.04);
  gStyle->SetPadTopMargin(0.06);
  ROOT::RDataFrame d0("DecayTree", kFull);
  auto d = d0.Define("ls", "D0_decayLength/D0_decayLengthErr")
               .Define("asym", "float p1=track_1_pT*cosh(track_1_pseudorapidity), p2=track_2_pT*cosh(track_2_pseudorapidity); return (p1-p2)/(p1+p2);")
               .Define("nV", "(int) nTracksOfVertex")
               .Define("nBC", "(int) nTracksOfBC")
               .Define("nPV", "(int) nPrimaryVerticesOfBC")
               .Filter("D0_decayLength<0.1");
  const std::string SIG = "D0_DIRA>0.995 && ls>2", CTL = "D0_DIRA>0.95 && D0_DIRA<0.98";
  auto ds = d.Filter(SIG), dc = d.Filter(CTL);

  // distributions and the PV-pull proxies vs multiplicity
  auto hnS = ds.Histo1D({"hnS", "", 60, 0, 120}, "nV"), hnC = dc.Histo1D({"hnC", "", 60, 0, 120}, "nV");
  auto hbS = ds.Histo1D({"hbS", "", 60, 0, 300}, "nBC");
  auto hpv = ds.Histo1D({"hpv", "", 10, 0, 10}, "nPV");
  auto pDca = d.Define("md", "std::min(abs(track_1_PV_DCA_sig),abs(track_2_PV_DCA_sig))").Profile1D({"pDca", "", 30, 0, 90}, "nV", "md");
  auto pLs = d.Profile1D({"pLs", "", 30, 0, 90}, "nV", "ls");
  auto pPVv = d.Define("pvr", "cbrt(primary_vertex_volume)*1e4").Profile1D({"pPVv", "", 30, 0, 90}, "nV", "pvr");
  auto nS = ds.Count();
  auto nEq = ds.Filter("nV==nBC").Count();
  auto nOne = ds.Filter("nPV==1").Count();

  struct Bin
  {
    std::string key, label, cut;
  };
  std::vector<Bin> bins = {
      {"All", "all", "1"},
      {"VA", "N_{vtx} < 10", "nV<10"},
      {"VB", "10 #leq N_{vtx} < 20", "nV>=10 && nV<20"},
      {"VC", "20 #leq N_{vtx} < 35", "nV>=20 && nV<35"},
      {"VD", "N_{vtx} #geq 35", "nV>=35"},
      {"BA", "N_{BC} < 20", "nBC<20"},
      {"BB", "20 #leq N_{BC} < 40", "nBC>=20 && nBC<40"},
      {"BC", "N_{BC} #geq 40", "nBC>=40"},
      {"PV1", "one PV in the crossing", "nPV==1"},
  };
  const double ptb[] = {1.5, 2.0, 2.5, 3.0, 4.0, 20.0};
  struct Book
  {
    ROOT::RDF::RResultPtr<TH2D> mS, mC;
    ROOT::RDF::RResultPtr<TH1D> hS;
    ROOT::RDF::RResultPtr<ULong64_t> nS, nC;
  };
  std::vector<Book> bk;
  for (auto &b : bins)
  {
    auto s = ds.Filter(b.cut), c = dc.Filter(b.cut);
    TH2D model(("m" + b.key).c_str(), "", 20, -1, 1, 5, ptb);
    bk.push_back({s.Histo2D(ROOT::RDF::TH2DModel(model), "asym", "D0_pT"), c.Histo2D(ROOT::RDF::TH2DModel(model), "asym", "D0_pT"),
                  s.Histo1D({("hS" + b.key).c_str(), "", nb, mlo, mhi}, "D0_mass"), s.Count(), c.Count()});
  }
  L("pointing (L<1 mm, %s): %llu; nTracksOfVertex == nTracksOfBC: %.1f%%; one PV in the crossing: %.1f%%\n", SIG.c_str(), *nS, 100. * *nEq / *nS,
    100. * *nOne / *nS);
  L("pointing: N_vtx median-ish mean %.1f (control %.1f); N_BC mean %.1f\n", hnS->GetMean(), hnC->GetMean(), hbS->GetMean());
  V("hfourNs", Form("%llu", *nS));
  V("hfourEq", Form("%.1f", 100. * *nEq / *nS));
  V("hfourOnePV", Form("%.1f", 100. * *nOne / *nS));
  V("hfourMeanNv", Form("%.1f", hnS->GetMean()));
  V("hfourMeanNvC", Form("%.1f", hnC->GetMean()));
  L("proxies vs N_vtx (all L<1 mm candidates): <min daughter DCA significance>, <L/sigma_L>, <PV volume^(1/3) [um]>\n");
  for (int i = 1; i <= 30; i += 3)
    L("  N_vtx %3.0f: %.2f  %.2f  %.0f\n", pDca->GetBinCenter(i), pDca->GetBinContent(i), pLs->GetBinContent(i), pPVv->GetBinContent(i));

  std::vector<double> xs, ys, eys;
  L("D0-window excess with the (A,pT)-weighted control, by multiplicity:\n");
  for (size_t i = 0; i < bins.size(); ++i)
  {
    auto *W = (TH2D *) bk[i].mS->Clone(Form("W%zu", i));
    W->Divide(bk[i].mC.GetPtr());
    auto look = [W](float a, float p) { return W->GetBinContent(W->FindFixBin(std::max(-0.9999f, std::min(a, 0.9999f)), std::min(p, 19.9f))); };
    auto P = dc.Filter(bins[i].cut).Define("w", look, {"asym", "D0_pT"}).Histo1D({Form("P%zu", i), "", nb, mlo, mhi}, "D0_mass", "w");
    Res r = analyse(bk[i].hS.GetPtr(), P.GetPtr());
    double sb = r.ex / std::sqrt(std::max(1.0, r.eex * r.eex - r.ex));
    L("  %-4s %-28s Nsig %8llu Nctl %8llu | excess %6.0f +- %4.0f (S/sqrt(B) %.2f) | narrow z %5.2f | broad dchi2 %5.1f\n", bins[i].key.c_str(),
      bins[i].label.c_str(), *bk[i].nS, *bk[i].nC, r.ex, r.eex, sb, r.z, r.dchi);
    V("hfour" + bins[i].key + "Ns", Form("%llu", *bk[i].nS));
    V("hfour" + bins[i].key + "Ex", Form("%.0f", r.ex));
    V("hfour" + bins[i].key + "eEx", Form("%.0f", r.eex));
    V("hfour" + bins[i].key + "Z", Form("%.1f", r.z));
    V("hfour" + bins[i].key + "Dchi", Form("%.0f", r.dchi));
    if (bins[i].key[0] == 'V')
    {
      xs.push_back(i);
      ys.push_back(r.ex);
      eys.push_back(r.eex);
    }
  }

  TCanvas c("c", "c", 1500, 520);
  c.Divide(3, 1);
  c.cd(1);
  gPad->SetLogy();
  hnC->Scale(1. / hnC->Integral());
  hnS->Scale(1. / hnS->Integral());
  hnS->SetLineColor(kRed + 1);
  hnS->SetLineWidth(2);
  hnC->SetLineWidth(2);
  hnC->GetXaxis()->SetTitle("N_{vtx} (tracks of the candidate's PV)");
  hnC->GetYaxis()->SetTitle("fraction");
  hnC->DrawCopy("hist");
  hnS->DrawCopy("hist same");
  TLegend lg(0.5, 0.78, 0.95, 0.92);
  lg.SetBorderSize(0);
  lg.SetTextSize(0.045);
  lg.AddEntry(hnS.GetPtr(), "pointing", "l");
  lg.AddEntry(hnC.GetPtr(), "non-pointing", "l");
  lg.DrawClone();
  c.cd(2);
  pDca->SetMarkerStyle(20);
  pDca->SetMarkerSize(0.7);
  pDca->GetXaxis()->SetTitle("N_{vtx}");
  pDca->GetYaxis()->SetTitle("#LT min daughter DCA_{PV}/#sigma #GT");
  pDca->DrawCopy("e");
  c.cd(3);
  auto *h = new TH1D("hex", "", 4, 0, 4);
  const char *lab[4] = {"<10", "10-20", "20-35", "#geq35"};
  for (int k = 0; k < 4; ++k)
  {
    h->SetBinContent(k + 1, ys[k]);
    h->SetBinError(k + 1, eys[k]);
    h->GetXaxis()->SetBinLabel(k + 1, lab[k]);
  }
  h->GetXaxis()->SetTitle("N_{vtx}");
  h->GetXaxis()->SetLabelSize(0.07);
  h->GetYaxis()->SetTitle("excess in 1.845-1.885 GeV/c^{2}");
  h->SetMarkerStyle(20);
  h->Draw("e");
  TLine l0(0, 0, 4, 0);
  l0.SetLineStyle(2);
  l0.Draw();
  c.SaveAs("figures/hump4_mult.pdf");
  c.SaveAs("figures/hump4_mult.png");
  fclose(gL);
  fclose(gV);
}
