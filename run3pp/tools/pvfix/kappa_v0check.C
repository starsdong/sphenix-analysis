// kappa_v0check.C : validate the momentum-scale map on K0S -> pi pi and Lambda -> p pi (prod10pct V0 ntuples).
// Mass recomputed from the stored daughter (pT, eta, phi) with the stored hypotheses, without (rc) and with (kc) the
// kappa correction; Gaussian + pol2 fits of the peak, in pT bins of the V0 and for the training runs (79509-79516)
// vs the other runs.
//   root -l -b -q 'kappa_v0check.C+("<map csv>", "<Ks file glob>", "<Lambda file glob>", "<out prefix>")'
#include "kappa_apply.C"

#include <TCanvas.h>
#include <TF1.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TStyle.h>

namespace
{
  struct Fit
  {
    double mu, emu, sig, esig, n;
  };
  Fit fitPeak(TH1D *h, double m0, double s0, double lo, double hi)
  {
    // background first: pol2 to the sidebands (|m - m0| > 6 s0), then Gaussian + pol2 with the background as start
    static double sM0, sW;
    sM0 = m0;
    sW = 6 * s0;
    TF1 fb("fb", [](double *x, double *p) {
      if (std::fabs(x[0] - sM0) < sW) TF1::RejectPoint();
      return p[0] + p[1] * (x[0] - sM0) + p[2] * (x[0] - sM0) * (x[0] - sM0);
    }, lo, hi, 3);
    fb.SetParameters(h->GetBinContent(1), 0, 0);
    h->Fit(&fb, "QR0");
    TF1 f("f", "gaus(0)+[3]+[4]*(x-[6])+[5]*(x-[6])^2", lo, hi);
    double pk = h->GetBinContent(h->FindFixBin(m0)) - fb.Eval(m0);
    f.SetParameters(std::max(pk, 1.0), m0, s0, fb.GetParameter(0), fb.GetParameter(1), fb.GetParameter(2), m0);
    f.FixParameter(6, m0);
    f.SetParLimits(1, m0 - 3 * s0, m0 + 3 * s0);
    f.SetParLimits(2, 0.2 * s0, 5 * s0);
    h->Fit(&f, "QR0");
    return {f.GetParameter(1), f.GetParError(1), std::fabs(f.GetParameter(2)), f.GetParError(2), f.GetParameter(0) * std::fabs(f.GetParameter(2)) * std::sqrt(2 * M_PI) / h->GetBinWidth(1)};
  }
}  // namespace

void kappa_v0check(const char *csv, const char *ksIn, const char *laIn, const char *prefix)
{
  loadMap(csv);
  ROOT::EnableImplicitMT(6);
  gStyle->SetOptStat(0);
  FILE *out = fopen(Form("%s.txt", prefix), "w");
  struct Ch
  {
    const char *name, *mother, *in;
    double m0, s0, lo, hi;
  };
  for (Ch c : {Ch{"K0S", "K_S0", ksIn, 0.4976, 0.006, 0.44, 0.56}, Ch{"Lambda", "Lambda0", laIn, 1.11568, 0.0025, 1.085, 1.150}})
  {
    ROOT::RDataFrame d("DecayTree", c.in);
    auto e = d.Define("k1", [](char q, float p, float et, float f) { return kappaOf(q > 0 ? 1 : -1, p, et, f); }, {"track_1_charge", "track_1_pT", "track_1_pseudorapidity", "track_1_phi"})
                 .Define("k2", [](char q, float p, float et, float f) { return kappaOf(q > 0 ? 1 : -1, p, et, f); }, {"track_2_charge", "track_2_pT", "track_2_pseudorapidity", "track_2_phi"})
                 .Define("m1", [](int p) { return std::abs(p) == 2212 ? 0.938272f : 0.13957f; }, {"track_1_PDG_ID"})
                 .Define("m2", [](int p) { return std::abs(p) == 2212 ? 0.938272f : 0.13957f; }, {"track_2_PDG_ID"})
                 .Define("p1c", "k1*track_1_pT")
                 .Define("p2c", "k2*track_2_pT")
                 .Define("mrc", massOf, {"track_1_pT", "track_1_pseudorapidity", "track_1_phi", "m1", "track_2_pT", "track_2_pseudorapidity", "track_2_phi", "m2"})
                 .Define("mkc", massOf, {"p1c", "track_1_pseudorapidity", "track_1_phi", "m1", "p2c", "track_2_pseudorapidity", "track_2_phi", "m2"})
                 .Define("vpt", std::string(c.mother) + "_pT")
                 .Define("train", "runNumber>=79509 && runNumber<=79516");
    const double ptb[] = {0, 0.5, 1.0, 1.5, 2.0, 3.0, 5.0, 20};
    struct B
    {
      std::string sel, label;
      ROOT::RDF::RResultPtr<TH1D> r, k, s;
    };
    std::vector<B> bs;
    auto book = [&](const std::string &sel, const std::string &label) {
      auto x = e.Filter(sel);
      std::string t = std::to_string(bs.size());
      bs.push_back({sel, label, x.Histo1D({("r" + t + c.name).c_str(), "", 120, c.lo, c.hi}, "mrc"), x.Histo1D({("k" + t + c.name).c_str(), "", 120, c.lo, c.hi}, "mkc"),
                    x.Histo1D({("s" + t + c.name).c_str(), "", 120, c.lo, c.hi}, std::string(c.mother) + "_mass")});
    };
    book("1", "all");
    book("train", "training runs 79509-79516");
    book("!train", "other runs");
    for (int i = 0; i < 7; ++i) book(Form("vpt>=%g && vpt<%g", ptb[i], ptb[i + 1]), Form("pT %.1f-%.1f", ptb[i], ptb[i + 1]));
    auto n = d.Count();
    fprintf(out, "== %s (%s): %llu candidates; Gaussian+pol2, mu and sigma in MeV (PDG %.2f)\n", c.name, c.in, *n, 1e3 * c.m0);
    fprintf(out, "  %-28s | %-30s | %-30s | %-30s | sigma kc/rc\n", "selection", "stored KFP mass", "recomputed, uncorrected (rc)", "recomputed, corrected (kc)");
    for (auto &b : bs)
    {
      Fit s = fitPeak(b.s.GetPtr(), c.m0, c.s0, c.lo, c.hi), r = fitPeak(b.r.GetPtr(), c.m0, c.s0, c.lo, c.hi), k = fitPeak(b.k.GetPtr(), c.m0, c.s0, c.lo, c.hi);
      auto f = [](Fit x) { return std::string(Form("mu %7.2f+-%4.2f sig %5.2f+-%4.2f", 1e3 * x.mu, 1e3 * x.emu, 1e3 * x.sig, 1e3 * x.esig)); };
      fprintf(out, "  %-28s | %s | %s | %s | %.3f\n", b.label.c_str(), f(s).c_str(), f(r).c_str(), f(k).c_str(), k.sig / r.sig);
    }
    TCanvas cv("cv", "", 1200, 500);
    cv.Divide(2, 1);
    for (int p = 0; p < 2; ++p)
    {
      cv.cd(p + 1);
      auto &b = bs[p + 1];
      b.r->SetLineColor(kBlack);
      b.k->SetLineColor(kRed + 1);
      b.r->SetTitle(Form("%s, %s;m (GeV/c^{2});candidates", c.name, b.label.c_str()));
      b.r->SetMaximum(1.15 * std::max(b.r->GetMaximum(), b.k->GetMaximum()));
      b.r->DrawCopy("hist");
      b.k->DrawCopy("hist same");
      TLegend lg(0.6, 0.75, 0.95, 0.9);
      lg.AddEntry(b.r.GetPtr(), "uncorrected", "l");
      lg.AddEntry(b.k.GetPtr(), "kappa corrected", "l");
      lg.DrawClone();
    }
    cv.SaveAs(Form("%s_%s.png", prefix, c.name));
  }
  fclose(out);
  printf("%s", "");
  system(Form("cat %s.txt", prefix));
}
