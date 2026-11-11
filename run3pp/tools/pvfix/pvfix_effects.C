// pvfix_effects.C : effect of fix 1 and fix 2 within d0-pvfix D0 ntuples (original and *_PVrefit branches per candidate)
//   root -l -b -q 'pvfix_effects.C+("<file glob>", "<out prefix>")'  -> <prefix>.txt, <prefix>_values.txt, <prefix>_*.png/pdf
#include <ROOT/RDataFrame.hxx>
#include <TCanvas.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TStyle.h>

#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

namespace
{
  FILE *gL, *gV;
  template <typename... A>
  void L(const char *f, A... a)
  {
    printf(f, a...);
    fprintf(gL, f, a...);
  }
  void V(const std::string &k, double v, const char *fmt = "%.4g") { fprintf(gV, "%s %s\n", k.c_str(), Form(fmt, v)); }
  double q(TH1D &h, double p)
  {
    double x;
    h.GetQuantiles(1, &x, &p);
    return x;
  }
}  // namespace

void pvfix_effects(const char *files, const char *prefix)
{
  gL = fopen(Form("%s.txt", prefix), "w");
  gV = fopen(Form("%s_values.txt", prefix), "w");
  ROOT::EnableImplicitMT(4);
  gStyle->SetOptStat(0);
  ROOT::RDataFrame d0("DecayTree", files);
  auto d = d0.Define("ok", "primary_vertex_refit_status==1")
               .Define("okAll", "primary_vertex_refitAll_status==1")
               .Define("nIn", "(track_1_inPV==1)+(track_2_inPV==1)")
               .Define("dAx", "1e4*(primary_vertex_refitAll_x-primary_vertex_x)")
               .Define("dAy", "1e4*(primary_vertex_refitAll_y-primary_vertex_y)")
               .Define("dAz", "1e4*(primary_vertex_refitAll_z-primary_vertex_z)")
               .Define("dRxy", "1e4*hypot(primary_vertex_refit_x-primary_vertex_refitAll_x, primary_vertex_refit_y-primary_vertex_refitAll_y)")
               .Define("dRz", "1e4*(primary_vertex_refit_z-primary_vertex_refitAll_z)")
               .Define("dSxy", "1e4*hypot(primary_vertex_refit_x-primary_vertex_x, primary_vertex_refit_y-primary_vertex_y)")
               .Define("dca1o", "abs(track_1_PV_DCA)")
               .Define("dca2o", "abs(track_2_PV_DCA)")
               .Define("dca1r", "abs(track_1_PV_DCA_PVrefit)")
               .Define("dca2r", "abs(track_2_PV_DCA_PVrefit)")
               .Define("dca1a", "abs(track_1_PV_DCA_PVrefitAll)")
               .Define("dca2a", "abs(track_2_PV_DCA_PVrefitAll)")
               .Define("lso", "D0_decayLength/D0_decayLengthErr")
               .Define("lsa", "D0_decayLength_PVrefitAll/D0_decayLengthErr_PVrefitAll")
               .Define("lsr", "D0_decayLength_PVrefit/D0_decayLengthErr_PVrefit");
  auto n = d.Count();
  auto nOk = d.Filter("ok").Count(), nOkAll = d.Filter("okAll").Count();
  auto nIn0 = d.Filter("nIn==0").Count(), nIn1 = d.Filter("nIn==1").Count(), nIn2 = d.Filter("nIn==2").Count();
  auto nNoOther = d.Filter("nTracksOfVertex_noDaughters<2").Count();
  // fix 1: ID vs crossing
  auto nIdEqCross = d.Filter("primary_vertex_ID==primary_vertex_crossing").Count();
  // algorithm closure
  auto hAx = d.Filter("okAll").Histo1D({"hAx", "", 2000, -1000, 1000}, "dAx");
  auto hAz = d.Filter("okAll").Histo1D({"hAz", "", 2000, -1000, 1000}, "dAz");
  auto nAexact = d.Filter("okAll && abs(dAx)<1 && abs(dAy)<1 && abs(dAz)<1").Count();
  // fix 2 shift
  std::vector<ROOT::RDF::RResultPtr<TH1D>> hxy, hz;
  std::vector<ROOT::RDF::RResultPtr<ULong64_t>> nk;
  for (int k = 0; k <= 2; ++k)
  {
    auto f = d.Filter(Form("ok && okAll && nIn==%d", k));
    hxy.push_back(f.Histo1D({Form("hxy%d", k), ";|#Deltar_{xy}| refit - refitAll (#mum);candidates", 200, 0, 500}, "dRxy"));
    hz.push_back(f.Histo1D({Form("hz%d", k), ";#Deltaz refit - refitAll (#mum);candidates", 200, -500, 500}, "dRz"));
    nk.push_back(f.Count());
  }
  // quantities, original vs refit (status ok)
  auto f = d.Filter("ok && okAll");
  auto nF = f.Count();
  auto h1o = f.Histo1D({"h1o", ";|DCA_{PV}| (#mum);", 400, 0, 0.2}, "dca1o"), h1r = f.Histo1D({"h1r", "", 400, 0, 0.2}, "dca1r"),
       h1a = f.Histo1D({"h1a", "", 400, 0, 0.2}, "dca1a");
  auto hlo = f.Histo1D({"hlo", ";L/#sigma_{L};", 400, -10, 30}, "lso"), hlr = f.Histo1D({"hlr", "", 400, -10, 30}, "lsr"),
       hla = f.Histo1D({"hla", "", 400, -10, 30}, "lsa");
  auto hdo = f.Histo1D({"hdo", ";DIRA;", 400, 0.9, 1.0}, "D0_DIRA"), hdr = f.Histo1D({"hdr", "", 400, 0.9, 1.0}, "D0_DIRA_PVrefit"),
       hda = f.Histo1D({"hda", "", 400, 0.9, 1.0}, "D0_DIRA_PVrefitAll");
  // selections: skim (both |DCA|>60 um, DIRA>0.95) and pointing selection b of the d0cutscan note
  struct Sel
  {
    std::string key, o, r, a;
  };
  std::vector<Sel> sels = {
      {"Skim", "dca1o>0.006 && dca2o>0.006 && D0_DIRA>0.95", "dca1r>0.006 && dca2r>0.006 && D0_DIRA_PVrefit>0.95",
       "dca1a>0.006 && dca2a>0.006 && D0_DIRA_PVrefitAll>0.95"},
      {"Point", "dca1o>0.006 && dca2o>0.006 && D0_DIRA>0.995 && lso>2", "dca1r>0.006 && dca2r>0.006 && D0_DIRA_PVrefit>0.995 && lsr>2",
       "dca1a>0.006 && dca2a>0.006 && D0_DIRA_PVrefitAll>0.995 && lsa>2"},
      {"Tight", "dca1o>0.006 && dca2o>0.006 && D0_DIRA>0.999 && lso>5", "dca1r>0.006 && dca2r>0.006 && D0_DIRA_PVrefit>0.999 && lsr>5",
       "dca1a>0.006 && dca2a>0.006 && D0_DIRA_PVrefitAll>0.999 && lsa>5"},
  };
  std::vector<ROOT::RDF::RResultPtr<ULong64_t>> so, sr, sa, sb, sab;
  for (auto &s : sels)
  {
    so.push_back(f.Filter(s.o).Count());
    sr.push_back(f.Filter(s.r).Count());
    sa.push_back(f.Filter(s.a).Count());
    sb.push_back(f.Filter("(" + s.o + ") && (" + s.a + ")").Count());   // original and refitAll
    sab.push_back(f.Filter("(" + s.a + ") && (" + s.r + ")").Count());  // refitAll and refit
  }

  L("candidates %llu; refit status ok %llu (%.2f%%), refitAll ok %llu (%.2f%%); fewer than 2 other tracks in the PV %llu (%.2f%%)\n", *n, *nOk, 100. * *nOk / *n,
    *nOkAll, 100. * *nOkAll / *n, *nNoOther, 100. * *nNoOther / *n);
  V("n", *n, "%.0f");
  V("okFrac", 100. * *nOk / *n, "%.2f");
  V("okAllFrac", 100. * *nOkAll / *n, "%.2f");
  V("noOtherFrac", 100. * *nNoOther / *n, "%.2f");
  L("fix 1: primary_vertex_ID == primary_vertex_crossing for %llu (%.2f%%)\n", *nIdEqCross, 100. * *nIdEqCross / *n);
  V("idEqCross", 100. * *nIdEqCross / *n, "%.2f");
  L("daughters attached to the PV: none %.2f%%, one %.2f%%, both %.2f%%\n", 100. * *nIn0 / *n, 100. * *nIn1 / *n, 100. * *nIn2 / *n);
  V("inNone", 100. * *nIn0 / *n, "%.2f");
  V("inOne", 100. * *nIn1 / *n, "%.2f");
  V("inBoth", 100. * *nIn2 / *n, "%.2f");
  L("algorithm closure refitAll - stored PV: within 1 um in x,y,z for %.2f%% of refitAll-ok; x rms %.1f um, z rms %.1f um, |dx| 99%% quantile %.1f um\n",
    100. * *nAexact / *nOkAll, hAx->GetRMS(), hAz->GetRMS(), [&] { TH1D h(*hAx); for (int i = 1; i <= h.GetNbinsX(); ++i) h.SetBinContent(i, hAx->GetBinContent(i) + hAx->GetBinContent(h.GetNbinsX() + 1 - i)); return q(h, 0.995); }() );
  V("closureExact", 100. * *nAexact / *nOkAll, "%.2f");
  for (int k = 0; k <= 2; ++k)
  {
    L("fix 2 shift, %d daughter(s) in PV (%llu cand.): |dr_xy| median %.1f um, 90%% %.1f um; dz rms %.1f um\n", k, *nk[k], q(*hxy[k], 0.5), q(*hxy[k], 0.9), hz[k]->GetRMS());
    V(Form("shiftXYmed%d", k), q(*hxy[k], 0.5), "%.1f");
    V(Form("shiftXY90%d", k), q(*hxy[k], 0.9), "%.1f");
    V(Form("shiftZrms%d", k), hz[k]->GetRMS(), "%.1f");
  }
  L("candidates with both refits: %llu (%.2f%%)\n", *nF, 100. * *nF / *n);
  V("nBoth", *nF, "%.0f");
  L("medians original -> refitAll -> refit: |DCA_1| %.1f -> %.1f -> %.1f um; L/sigma %.2f -> %.2f -> %.2f; DIRA %.5f -> %.5f -> %.5f\n", 1e4 * q(*h1o, 0.5),
    1e4 * q(*h1a, 0.5), 1e4 * q(*h1r, 0.5), q(*hlo, 0.5), q(*hla, 0.5), q(*hlr, 0.5), q(*hdo, 0.5), q(*hda, 0.5), q(*hdr, 0.5));
  V("dcaMedO", 1e4 * q(*h1o, 0.5), "%.1f");
  V("dcaMedA", 1e4 * q(*h1a, 0.5), "%.1f");
  V("dcaMedR", 1e4 * q(*h1r, 0.5), "%.1f");
  V("lsMedO", q(*hlo, 0.5), "%.2f");
  V("lsMedA", q(*hla, 0.5), "%.2f");
  V("lsMedR", q(*hlr, 0.5), "%.2f");
  for (size_t i = 0; i < sels.size(); ++i)
  {
    L("selection %-6s original %llu | refitAll %llu (x%.3f; kept %.1f%% of original) | refit %llu (x%.3f vs refitAll; kept %.1f%% of refitAll, new %.1f%%)\n",
      sels[i].key.c_str(), *so[i], *sa[i], (double) *sa[i] / *so[i], 100. * *sb[i] / *so[i], *sr[i], (double) *sr[i] / *sa[i], 100. * *sab[i] / *sa[i],
      100. * (*sr[i] - *sab[i]) / *sa[i]);
    V("sel" + sels[i].key + "O", *so[i], "%.0f");
    V("sel" + sels[i].key + "A", *sa[i], "%.0f");
    V("sel" + sels[i].key + "R", *sr[i], "%.0f");
    V("sel" + sels[i].key + "AO", (double) *sa[i] / *so[i], "%.3f");
    V("sel" + sels[i].key + "RA", (double) *sr[i] / *sa[i], "%.3f");
  }
  // figures
  TCanvas c("c", "c", 1500, 500);
  c.Divide(3, 1);
  auto pair = [&](int pad, TH1D *o, TH1D *r, TH1D *a, bool logy) {
    c.cd(pad);
    if (logy) gPad->SetLogy();
    o->SetLineColor(kBlack);
    r->SetLineColor(kRed + 1);
    o->SetLineWidth(2);
    r->SetLineWidth(2);
    o->SetMaximum(std::max(o->GetMaximum(), r->GetMaximum()) * (logy ? 3 : 1.2));
    a->SetLineColor(kBlue + 1);
    a->SetLineWidth(2);
    a->SetLineStyle(2);
    o->DrawCopy("hist");
    a->DrawCopy("hist same");
    r->DrawCopy("hist same");
    TLegend lg(0.45, 0.74, 0.9, 0.9);
    lg.SetBorderSize(0);
    lg.AddEntry(o, "original PV", "l");
    lg.AddEntry(a, "refit, all tracks", "l");
    lg.AddEntry(r, "refit w/o daughters", "l");
    lg.DrawClone();
  };
  h1o->GetXaxis()->SetTitle("|DCA_{PV}| track 1 (cm)");
  pair(1, h1o.GetPtr(), h1r.GetPtr(), h1a.GetPtr(), true);
  pair(2, hlo.GetPtr(), hlr.GetPtr(), hla.GetPtr(), true);
  pair(3, hdo.GetPtr(), hdr.GetPtr(), hda.GetPtr(), true);
  c.SaveAs(Form("%s_quantities.png", prefix));
  c.SaveAs(Form("%s_quantities.pdf", prefix));
  TCanvas c2("c2", "c2", 1000, 450);
  c2.Divide(2, 1);
  int col[3] = {kBlack, kBlue + 1, kRed + 1};
  TLegend lg(0.5, 0.7, 0.9, 0.9);
  lg.SetBorderSize(0);
  for (int k = 0; k <= 2; ++k)
  {
    for (int p = 0; p < 2; ++p)
    {
      c2.cd(p + 1);
      gPad->SetLogy();
      TH1D *h = p == 0 ? hxy[k].GetPtr() : hz[k].GetPtr();
      h->SetLineColor(col[k]);
      h->SetLineWidth(2);
      if (h->Integral() > 0) h->Scale(1. / h->Integral());
      h->SetMaximum(1);
      h->SetMinimum(1e-5);
      h->DrawCopy(k ? "hist same" : "hist");
    }
    lg.AddEntry(hxy[k].GetPtr(), Form("%d daughter(s) in PV", k), "l");
  }
  c2.cd(1);
  lg.DrawClone();
  c2.SaveAs(Form("%s_shift.png", prefix));
  c2.SaveAs(Form("%s_shift.pdf", prefix));
  fclose(gL);
  fclose(gV);
}
