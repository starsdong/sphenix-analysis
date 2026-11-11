// pvfix_d0_v2.C : D0 signal search with three PV versions (original [stock in-job selection], pair refit w/o daughters,
// KFParticle fit w/o daughters), 10% d0-pvfix v2 skim. Derived from pvfix_d0.C.
//   root -l -b -q 'pvfix_d0_v2.C+("<skim file>", "<out prefix>")'  -> <prefix>.txt, <prefix>_values.txt, <prefix>_*.png/pdf
// Per version v and pT range: prompt-like (L_v < 1 mm), pointing = DIRA_v>0.995 && L_v/sigma>2, tight = DIRA_v>0.999 && L_v/sigma>5,
// control = DIRA_v 0.95-0.98; both daughters |DCA_v| > 60 um. Background model = control weighted with the pointing probability
// measured in bins of (momentum asymmetry A, D0 pT) (d0cutscan note, Sec. 9). Window excess in 1.845-1.885 GeV/c^2 above a pol2
// fitted to data/prediction outside the window; injection test with Ninj control candidates from 1.855-1.875.
// Base sample: candidates for which both refits w/o daughters exist; the original version additionally requires
// pvfix_passOriginal (in-job stock cuts), the refit versions do not (that is the point of v2).
#include <ROOT/RDataFrame.hxx>
#include <TCanvas.h>
#include <TF1.h>
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
    if (nrm && r->Integral() > 0) r->Scale(r->GetNbinsX() / r->Integral());
    return r;
  }
  struct Res
  {
    double ex, eex, z, dchi, chi2;
    int ndf;
  };
  Res analyse(const TH1D *S, const TH1D *P, const char *tag)
  {
    TH1D *dd = hratio(S, P, Form("dd_%s", tag), false);
    double k = dd->GetNbinsX() / dd->Integral();
    dd->Scale(k);
    TF1 f(Form("fx_%s", tag), [](double *x, double *p) {
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
    return {ex, std::sqrt(n), g.GetParError(0) > 0 ? g.GetParameter(0) / g.GetParError(0) : 0, p2.GetChisquare() - b.GetChisquare(), p2.GetChisquare(),
            p2.GetNDF()};
  }
}  // namespace

void pvfix_d0_v2(const char *in, const char *prefix, const char *extra = "1", const char *extraLabel = "")
{
  gL = fopen(Form("%s.txt", prefix), "w");
  gV = fopen(Form("%s_values.txt", prefix), "w");
  ROOT::EnableImplicitMT(6);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetTitleSize(0.05, "XYZ");
  gStyle->SetLabelSize(0.045, "XYZ");
  gStyle->SetPadLeftMargin(0.16);
  gStyle->SetPadBottomMargin(0.13);
  gStyle->SetPadRightMargin(0.04);
  gStyle->SetPadTopMargin(0.07);
  ROOT::RDataFrame d0("DecayTree", in);
  auto base = d0.Filter("primary_vertex_refit_status==1 && primary_vertex_kf_status==1")
                  .Define("dz1", "float a=track_1_PV_DCA, b=track_1_PV_DCA_xy; return (float) std::sqrt(std::max(0.f, a*a-b*b));")
                  .Define("dz2", "float a=track_2_PV_DCA, b=track_2_PV_DCA_xy; return (float) std::sqrt(std::max(0.f, a*a-b*b));")
                  .Define("asym", "float p1=track_1_pT*cosh(track_1_pseudorapidity), p2=track_2_pT*cosh(track_2_pseudorapidity); return (p1-p2)/(p1+p2);");
  auto nAll = d0.Count(), nBase = base.Count();
  struct Ver
  {
    std::string key, sfx, label;
  };
  std::vector<Ver> vers = {{"O", "", "original PV"}, {"R", "_PVrefit", "pair refit w/o daughters"}, {"K", "_PVkf", "KF fit w/o daughters"}};
  struct Rng
  {
    std::string key, cut, label;
  };
  std::vector<Rng> rngs = {{"Hi", "D0_pT>1.5", "p_{T} > 1.5"}, {"All", "1", "all p_{T}"}};
  struct Sel
  {
    std::string key, label;
  };
  std::vector<Sel> sels = {{"P", "DIRA>0.995, L/#sigma>2"}, {"T", "DIRA>0.999, L/#sigma>5"}};
  const double ptb[] = {0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 20.0};
  struct Book
  {
    std::string tag, vk, rk, sk;
    ROOT::RDF::RResultPtr<TH2D> mS, mC, mI;
    ROOT::RDF::RResultPtr<TH1D> hS, hC, hI;
    ROOT::RDF::RNode dc;
  };
  std::vector<Book> bk;
  for (auto &v : vers)
    for (auto &r : rngs)
      for (auto &s : sels)
      {
        std::string dca = "abs(track_1_PV_DCA" + v.sfx + ")>0.006 && abs(track_2_PV_DCA" + v.sfx + ")>0.006";
        std::string ls = "(D0_decayLength" + v.sfx + "/D0_decayLengthErr" + v.sfx + ")";
        std::string prompt = "D0_decayLength" + v.sfx + "<0.1";
        std::string sig = s.key == "P" ? "D0_DIRA" + v.sfx + ">0.995 && " + ls + ">2" : "D0_DIRA" + v.sfx + ">0.999 && " + ls + ">5";
        std::string ctl = "D0_DIRA" + v.sfx + ">0.95 && D0_DIRA" + v.sfx + "<0.98";
        std::string common = std::string("(") + extra + ") && " + r.cut + " && " + dca + " && " + prompt + (v.key == "O" ? std::string(" && pvfix_passOriginal==1") : std::string(""));
        auto ds = base.Filter(common + " && " + sig);
        auto dc = base.Filter(common + " && " + ctl);
        auto inj = dc.Filter("D0_mass>1.855 && D0_mass<1.875");
        std::string tag = v.key + r.key + s.key;
        TH2D model(("m" + tag).c_str(), "", 20, -1, 1, 7, ptb);
        bk.push_back({tag, v.key, r.key, s.key, ds.Histo2D(ROOT::RDF::TH2DModel(model), "asym", "D0_pT"), dc.Histo2D(ROOT::RDF::TH2DModel(model), "asym", "D0_pT"),
                      inj.Histo2D(ROOT::RDF::TH2DModel(model), "asym", "D0_pT"), ds.Histo1D({("hS" + tag).c_str(), "", nb, mlo, mhi}, "D0_mass"),
                      dc.Histo1D({("hC" + tag).c_str(), "", nb, mlo, mhi}, "D0_mass"), inj.Histo1D({("hI" + tag).c_str(), "", nb, mlo, mhi}, "D0_mass"), dc});
      }
  L("extra selection: %s [%s]\n", extra, extraLabel);
  auto nExtra = base.Filter(extra).Count();
  L("skim candidates %llu, with pair refit and KF fit %llu (%.2f%%), passing the extra selection %llu (%.2f%% of these)\n", *nAll, *nBase, 100. * *nBase / *nAll, *nExtra, 100. * *nExtra / *nBase);
  V("skimN", Form("%llu", *nAll));
  V("baseN", Form("%llu", *nBase));
  std::vector<TH1D *> rawR, dblR;
  for (auto &b : bk)
  {
    double nS = b.hS->GetEntries(), nC = b.hC->GetEntries();
    const double Ninj = std::max(50.0, std::round(0.002 * nS));  // 0.2% of the pointing sample
    auto *W = (TH2D *) b.mS->Clone(("W" + b.tag).c_str());
    W->Divide(b.mC.GetPtr());
    auto *Si = (TH2D *) b.mS->Clone(("Si" + b.tag).c_str());
    Si->Add(b.mI.GetPtr(), Ninj / std::max(1.0, b.mI->GetEntries()));
    auto *Wi = (TH2D *) Si->Clone(("Wi" + b.tag).c_str());
    Wi->Divide(b.mC.GetPtr());
    auto look = [](const TH2D *w) { return [w](float a, float p) { return w->GetBinContent(w->FindFixBin(std::max(-0.9999f, std::min(a, 0.9999f)), std::min(p, 19.9f))); }; };
    auto P = b.dc.Define("w", look(W), {"asym", "D0_pT"}).Histo1D({("P" + b.tag).c_str(), "", nb, mlo, mhi}, "D0_mass", "w");
    auto Pi = b.dc.Define("wi", look(Wi), {"asym", "D0_pT"}).Histo1D({("Pi" + b.tag).c_str(), "", nb, mlo, mhi}, "D0_mass", "wi");
    auto *hSi = (TH1D *) b.hS->Clone(("hSi" + b.tag).c_str());
    hSi->Add(b.hI.GetPtr(), Ninj / std::max(1.0, b.hI->GetEntries()));
    TH1D *R = hratio(b.hS.GetPtr(), b.hC.GetPtr(), ("R" + b.tag).c_str());
    TF1 p2("p2r", "pol2", 1.65, 2.08);
    R->Fit(&p2, "QR0");
    TF1 bg("bgr", "[0]*exp(-0.5*((x-[1])/[2])^2)+pol2(3)", 1.65, 2.08);
    bg.SetParameters(0.05, 1.88, 0.08, p2.GetParameter(0), p2.GetParameter(1), p2.GetParameter(2));
    bg.SetParLimits(1, 1.75, 2.0);
    bg.SetParLimits(2, 0.005, 0.25);
    R->Fit(&bg, "QR0");
    double rawDchi = p2.GetChisquare() - bg.GetChisquare();
    Res r0 = analyse(b.hS.GetPtr(), P.GetPtr(), b.tag.c_str()), r1 = analyse(hSi, Pi.GetPtr(), (b.tag + "i").c_str());
    double rec = r1.ex - r0.ex;
    L("%-4s %-22s %-11s %-24s Nsig %8.0f Nctl %8.0f | raw ratio broad dchi2 %5.0f | weighted: pol2 chi2/ndf %4.0f/%d, broad dchi2 %4.0f | excess %6.0f +- %4.0f "
      "(%.1f sigma), narrow z %5.2f | injected %.0f -> recovered %.0f%%\n",
      b.tag.c_str(), vers[b.vk == "O" ? 0 : b.vk == "R" ? 1 : 2].label.c_str(), b.rk == "Hi" ? "pT>1.5" : "all pT", b.sk == "P" ? "DIRA>0.995,L/s>2" : "DIRA>0.999,L/s>5",
      nS, nC, rawDchi, r0.chi2, r0.ndf, r0.dchi, r0.ex, r0.eex, r0.ex / r0.eex, r0.z, Ninj, 100 * rec / Ninj);
    V("d" + b.tag + "Ns", Form("%.0f", nS));
    V("d" + b.tag + "Nc", Form("%.0f", nC));
    V("d" + b.tag + "RawDchi", Form("%.0f", rawDchi));
    V("d" + b.tag + "Dchi", Form("%.0f", r0.dchi));
    V("d" + b.tag + "Ex", Form("%.0f", r0.ex));
    V("d" + b.tag + "eEx", Form("%.0f", r0.eex));
    V("d" + b.tag + "Z", Form("%.1f", r0.z));
    V("d" + b.tag + "Rec", Form("%.0f", 100 * rec / Ninj));
    rawR.push_back(R);
    dblR.push_back(hratio(b.hS.GetPtr(), P.GetPtr(), ("D" + b.tag).c_str()));
  }
  // figures: for each (pT range, selection): mass spectra, raw ratio, double ratio; three versions overlaid
  int col[3] = {kBlack, kBlue + 1, kRed + 1};
  for (auto &r : rngs)
    for (auto &s : sels)
    {
      TCanvas c("c", "c", 1500, 500);
      c.Divide(3, 1);
      TLegend lg(0.2, 0.72, 0.7, 0.92);
      lg.SetBorderSize(0);
      lg.SetTextSize(0.042);
      for (int iv = 0; iv < 3; ++iv)
      {
        size_t i = 0;
        for (; i < bk.size(); ++i)
          if (bk[i].vk == vers[iv].key && bk[i].rk == r.key && bk[i].sk == s.key) break;
        auto *h = (TH1D *) bk[i].hS->Clone(Form("hm%zu", i));
        h->SetLineColor(col[iv]);
        h->SetMarkerColor(col[iv]);
        h->SetMarkerStyle(20 + iv);
        h->SetMarkerSize(0.6);
        h->GetXaxis()->SetTitle("m(K#pi) (GeV/c^{2})");
        h->GetYaxis()->SetTitle("candidates / 10 MeV");
        c.cd(1);
        h->SetMinimum(0);
        h->DrawCopy(iv ? "e same" : "e");
        lg.AddEntry(h, vers[iv].label.c_str(), "pe");
        for (int p = 0; p < 2; ++p)
        {
          TH1D *x = p == 0 ? rawR[i] : dblR[i];
          x->SetLineColor(col[iv]);
          x->SetMarkerColor(col[iv]);
          x->SetMarkerStyle(20 + iv);
          x->SetMarkerSize(0.6);
          x->GetXaxis()->SetTitle("m(K#pi) (GeV/c^{2})");
          x->GetYaxis()->SetTitle(p == 0 ? "pointing / non-pointing (norm.)" : "data / weighted prediction (norm.)");
          x->SetMinimum(p == 0 ? 0.6 : 0.85);
          x->SetMaximum(p == 0 ? 1.5 : 1.15);
          c.cd(p + 2);
          x->DrawCopy(iv ? "e same" : "e");
        }
      }
      for (int p = 1; p <= 3; ++p)
      {
        c.cd(p);
        gPad->Update();
        TLine l(mD0, gPad->GetUymin(), mD0, gPad->GetUymax());
        l.SetLineStyle(3);
        l.DrawClone();
      }
      c.cd(1);
      lg.DrawClone();
      TLatex t;
      t.SetNDC();
      t.SetTextSize(0.045);
      t.DrawLatex(0.2, 0.95, Form("%s, %s, L < 1 mm, 10%% sample %s", r.label.c_str(), s.label.c_str(), extraLabel));
      c.SaveAs(Form("%s_%s%s.png", prefix, r.key.c_str(), s.key.c_str()));
      c.SaveAs(Form("%s_%s%s.pdf", prefix, r.key.c_str(), s.key.c_str()));
    }
  fclose(gL);
  fclose(gV);
}
