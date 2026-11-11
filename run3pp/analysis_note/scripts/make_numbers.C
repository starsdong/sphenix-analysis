// make_numbers.C : all numbers, tables and figures for the run3pp KFParticle analysis note.
//
// usage (sPHENIX env, ana.581):
//   root -b -q 'make_numbers.C+("<merged ntuple dir>", "<note root>")'
// reads   <dir>/KFP_{EventInfo,Ks,Lambda,D0,Ks_noXmatch,Lambda_noXmatch}_run79509_seg000-030.root
// writes  <note>/logs/numbers.txt, <note>/tables/*.tex, <note>/figures/*.{pdf,png}
#include <TCanvas.h>
#include <TF1.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TTreeReaderValue.h>

#include <cmath>
#include <cstdio>
#include <fstream>
#include <map>
#include <set>
#include <string>
#include <tuple>
#include <vector>

namespace
{
  std::string gNote, gDir, gTag = "run79509_seg000-030";
  FILE *gLog = nullptr;
  std::ofstream gMac;  // tables/numbers.tex

  // macro names: letters only
  void num(const std::string &name, const std::string &val)
  {
    fprintf(gLog, "%s = %s\n", name.c_str(), val.c_str());
    gMac << "\\newcommand{\\" << name << "}{" << val << "}\n";
  }
  std::string F(double v, int p = 1)
  {
    char b[64];
    snprintf(b, 64, "%.*f", p, v);
    return b;
  }
  std::string I(double v)
  {
    // thousands separator
    char b[64];
    snprintf(b, 64, "%.0f", v);
    std::string s = b, o;
    int n = 0;
    for (int i = (int) s.size() - 1; i >= 0; --i)
    {
      o.insert(o.begin(), s[i]);
      if (++n % 3 == 0 && i > 0 && s[i - 1] != '-') o.insert(o.begin(), ',');
    }
    return o;
  }
  void save(TCanvas &c, const std::string &name)
  {
    c.SaveAs((gNote + "/figures/" + name + ".pdf").c_str());
    c.SaveAs((gNote + "/figures/" + name + ".png").c_str());
  }
  TLatex lat()
  {
    TLatex l;
    l.SetNDC();
    l.SetTextSize(0.045);
    return l;
  }
  void style()
  {
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetPadLeftMargin(0.17);
    gStyle->SetPadBottomMargin(0.14);
    gStyle->SetPadRightMargin(0.04);
    gStyle->SetPadTopMargin(0.05);
    gStyle->SetLabelSize(0.045, "XY");
    gStyle->SetTitleSize(0.05, "XY");
    gStyle->SetTitleOffset(1.6, "Y");
    gStyle->SetLegendBorderSize(0);
    gStyle->SetEndErrorSize(0);
  }

  char dtype(unsigned nintt, unsigned nmvtx) { return nintt > 0 ? 'I' : (nmvtx > 0 ? 'M' : 'T'); }
  std::string pairType(char a, char b)
  {
    std::string order = "IMT";
    if (order.find(a) > order.find(b)) std::swap(a, b);
    return std::string(1, a) + b;
  }

  struct Cand
  {
    float mass, dira, dl, dlerr, dlxy, dca, d1xy, d2xy, x, y, z, pt, chi2;
    unsigned ndof;
    int x1, x2, xpv, q;
    char t1, t2;
    Long64_t cbco;
    int evt;
    float p1x;
  };

  std::vector<Cand> readCands(const std::string &file, const std::string &M)
  {
    std::vector<Cand> v;
    TFile f(file.c_str());
    TTree *t = (TTree *) f.Get("DecayTree");
    TTreeReader r(t);
    TTreeReaderValue<float> mass(r, (M + "_mass").c_str()), dira(r, (M + "_DIRA").c_str()), dl(r, (M + "_decayLength").c_str()),
        dlerr(r, (M + "_decayLengthErr").c_str()), dlxy(r, (M + "_decayLength_xy").c_str()), x(r, (M + "_x").c_str()),
        y(r, (M + "_y").c_str()), z(r, (M + "_z").c_str()), pt(r, (M + "_pT").c_str()), chi2(r, (M + "_chi2").c_str()),
        dca(r, "track_1_track_2_DCA"), d1xy(r, "track_1_PV_DCA_xy"), d2xy(r, "track_2_PV_DCA_xy"), p1x(r, "track_1_px");
    TTreeReaderValue<unsigned> ndof(r, (M + "_nDoF").c_str());
    TTreeReaderValue<int> x1(r, "track_1_bunch_crossing"), x2(r, "track_2_bunch_crossing"), xpv(r, "primary_vertex_ID"),
        pdg1(r, "track_1_PDG_ID"), pdg2(r, "track_2_PDG_ID"), evt(r, "eventNumber");
    TTreeReaderValue<unsigned> i1(r, "track_1_INTT_nHits"), m1(r, "track_1_MVTX_nHits"), i2(r, "track_2_INTT_nHits"), m2(r, "track_2_MVTX_nHits");
    TTreeReaderValue<Long64_t> cbco(r, "Collision_BCO");
    while (r.Next())
    {
      Cand k;
      k.mass = *mass; k.dira = *dira; k.dl = *dl; k.dlerr = *dlerr; k.dlxy = *dlxy; k.dca = *dca;
      k.d1xy = std::fabs(*d1xy); k.d2xy = std::fabs(*d2xy); k.x = *x; k.y = *y; k.z = *z; k.pt = *pt;
      k.chi2 = *chi2; k.ndof = *ndof; k.x1 = *x1; k.x2 = *x2; k.xpv = *xpv;
      k.t1 = dtype(*i1, *m1); k.t2 = dtype(*i2, *m2);
      k.q = 0;
      if (std::abs(*pdg1) == 2212) k.q = (*pdg1 > 0 ? 1 : -1);
      if (std::abs(*pdg2) == 2212) k.q = (*pdg2 > 0 ? 1 : -1);
      k.cbco = *cbco; k.evt = *evt; k.p1x = *p1x;
      v.push_back(k);
    }
    return v;
  }

  // ---------------------------------------------------------------- QA (offline) selection for V0 peaks
  const double cDIRA = 0.995, cLxy = 0.5, cDCA = 0.2, cTrkDCAxy = 0.02;
  bool qa(const Cand &k) { return k.dira > cDIRA && k.dlxy > cLxy && k.dca < cDCA && k.d1xy > cTrkDCAxy && k.d2xy > cTrkDCAxy; }

  // duplicate flag: same collision BCO and kinematics already seen
  std::vector<bool> dupFlags(const std::vector<Cand> &v)
  {
    std::set<std::tuple<Long64_t, long, long>> seen;
    std::vector<bool> d(v.size(), false);
    for (size_t i = 0; i < v.size(); ++i)
    {
      auto key = std::make_tuple(v[i].cbco, std::lround(v[i].mass * 1e5), std::lround(v[i].p1x * 1e4));
      if (!seen.insert(key).second) d[i] = true;
    }
    return d;
  }

  struct Fit
  {
    double mean = 0, emean = 0, sigma = 0, esigma = 0, yield = 0, eyield = 0, bkg3 = 0;
    bool ok = false;
  };
  // Gaussian + pol2. If fixMean/fixSigma > 0 the shape is fixed and only the yield floats.
  Fit fitPeak(TH1F *h, double m0, double s0, double lo, double hi, bool fixShape, int color = kRed)
  {
    Fit r;
    if (h->Integral() < 20) return r;
    TF1 *f = new TF1(Form("f_%s", h->GetName()), "[0]*exp(-0.5*((x-[1])/[2])^2)+[3]+[4]*(x-[1])+[5]*(x-[1])^2", lo, hi);
    double bw = h->GetBinWidth(1);
    double b0 = 0.5 * (h->GetBinContent(h->FindBin(lo + bw)) + h->GetBinContent(h->FindBin(hi - bw)));
    double a0 = std::max(1.0, h->GetBinContent(h->FindBin(m0)) - b0);
    f->SetParameters(a0, m0, s0, b0, 0, 0);
    if (fixShape)
    {
      // amplitude unbounded (may go negative), mean and width fixed to the reference fit
      f->FixParameter(1, m0);
      f->FixParameter(2, s0);
    }
    else
    {
      f->SetParLimits(0, 0, 1e9);
      f->SetParLimits(1, m0 - 4 * s0, m0 + 4 * s0);
      f->SetParLimits(2, 0.3 * s0, 3 * s0);
    }
    f->SetLineColor(color);
    f->SetLineWidth(2);
    f->SetNpx(400);
    int st = h->Fit(f, "QRLS0");
    r.ok = (st == 0);
    if (!r.ok) fprintf(gLog, "# fit status %d for %s\n", st, h->GetName());
    r.mean = f->GetParameter(1);
    r.emean = fixShape ? 0 : f->GetParError(1);
    r.sigma = std::fabs(f->GetParameter(2));
    r.esigma = fixShape ? 0 : f->GetParError(2);
    r.yield = f->GetParameter(0) * r.sigma * std::sqrt(2 * M_PI) / bw;
    r.eyield = f->GetParError(0) * r.sigma * std::sqrt(2 * M_PI) / bw;
    TF1 p("p", "[0]+[1]*(x-[3])+[2]*(x-[3])^2", lo, hi);
    p.SetParameters(f->GetParameter(3), f->GetParameter(4), f->GetParameter(5), r.mean);
    r.bkg3 = p.Integral(r.mean - 3 * r.sigma, r.mean + 3 * r.sigma) / bw;
    h->GetListOfFunctions()->Add(f);
    return r;
  }

  void drawMass(TH1F *h, const std::string &xt, int color = kBlack)
  {
    h->SetMarkerStyle(20);
    h->SetMarkerSize(0.7);
    h->SetMarkerColor(color);
    h->SetLineColor(color);
    h->GetXaxis()->SetTitle(xt.c_str());
    h->GetYaxis()->SetTitle(Form("candidates / %.1f MeV", 1e3 * h->GetBinWidth(1)));
    h->SetMinimum(0);
    h->Draw("E");
    for (auto *o : *h->GetListOfFunctions())
      if (o->InheritsFrom("TF1")) ((TF1 *) o)->Draw("same");
  }
}  // namespace

void make_numbers(const std::string &dir, const std::string &note)
{
  style();
  TH1::SetDefaultSumw2();
  gDir = dir;
  gNote = note;
  gSystem->mkdir((note + "/logs").c_str(), true);
  gSystem->mkdir((note + "/tables").c_str(), true);
  gSystem->mkdir((note + "/figures").c_str(), true);
  gLog = fopen((note + "/logs/numbers.txt").c_str(), "w");
  gMac.open(note + "/tables/numbers.tex");
  gMac << "% generated by scripts/make_numbers.C -- do not edit\n";
  auto fn = [&](const std::string &ch) { return dir + "/KFP_" + ch + "_" + gTag + ".root"; };

  // ================================================================ events
  double nEv = 0;
  {
    TFile f(fn("EventInfo").c_str());
    TTree *t = (TTree *) f.Get("EventInfo");
    nEv = t->GetEntries();
    TTreeReader r(t);
    TTreeReaderValue<int> ntrk(r, "nTrk"), nI(r, "nTrkINTT"), nM(r, "nTrkMVTXonly"), nT(r, "nTrkTPConly"), nX0(r, "nTrkX0"), nV(r, "nVtx"), run(r, "runNumber");
    TTreeReaderArray<int> vX(r, "vtx_crossing");
    TTreeReaderArray<float> vz(r, "vtx_z");
    double s[6] = {0};
    long evX0 = 0, evX0multi = 0, nvtx = 0, nvtxX0 = 0, nNoVtx = 0;
    int xmin = 9999, xmax = -9999;
    std::set<int> runs;
    TH1F *hX = new TH1F("hVtxX", ";vertex beam crossing (relative to GL1 BCO);vertices", 600, -150.5, 449.5);
    TH1F *hZ = new TH1F("hVtxZ", ";vertex z (cm);vertices", 120, -30, 30);
    TH1F *hNtrk = new TH1F("hNtrk", ";tracks per event;events", 100, 0, 300);
    TH1F *hNtpc = new TH1F("hNtpc", ";TPC-only tracks per event;events", 100, 0, 100);
    while (r.Next())
    {
      s[0] += *ntrk; s[1] += *nI; s[2] += *nM; s[3] += *nT; s[4] += *nX0; s[5] += *nV;
      runs.insert(*run);
      hNtrk->Fill(*ntrk);
      hNtpc->Fill(*nT);
      if (*nV == 0) ++nNoVtx;
      int n0 = 0;
      for (size_t i = 0; i < vX.GetSize(); ++i)
      {
        ++nvtx;
        hX->Fill(vX[i]);
        hZ->Fill(vz[i]);
        xmin = std::min(xmin, vX[i]);
        xmax = std::max(xmax, vX[i]);
        if (vX[i] == 0) { ++n0; ++nvtxX0; }
      }
      if (n0 > 0) ++evX0;
      if (n0 > 1) ++evX0multi;
    }
    num("nEvents", I(nEv));
    num("nRuns", I(runs.size()));
    num("trkPerEv", F(s[0] / nEv));
    num("trkINTTPerEv", F(s[1] / nEv));
    num("trkMVTXonlyPerEv", F(s[2] / nEv));
    num("trkTPConlyPerEv", F(s[3] / nEv));
    num("trkTPConlyFrac", F(100 * s[3] / s[0]));
    num("trkXzeroPerEv", F(s[4] / nEv));
    num("vtxPerEv", F(s[5] / nEv));
    num("evNoVtxFrac", F(100. * nNoVtx / nEv));
    num("evVtxXzeroFrac", F(100. * evX0 / nEv));
    num("evVtxXzeroMultiFrac", F(100. * evX0multi / nEv));
    num("vtxXmin", I(xmin));
    num("vtxXmax", I(xmax));
    // per-crossing vertex density away from 0
    double other = 0;
    int nb = 0;
    for (int b = 1; b <= hX->GetNbinsX(); ++b)
    {
      double c = hX->GetBinCenter(b);
      if (c >= 20 && c <= 380) { other += hX->GetBinContent(b); ++nb; }
    }
    num("vtxAtXzero", I(hX->GetBinContent(hX->FindBin(0))));
    num("vtxPerOtherX", F(other / nb));

    TCanvas c("c", "", 1500, 500);
    c.Divide(3, 1);
    c.cd(1);
    gPad->SetLogy();
    hX->SetLineColor(kBlue + 1);
    hX->Draw("hist");
    auto l = lat();
    l.DrawLatex(0.5, 0.88, Form("%s events", I(nEv).c_str()));
    c.cd(2);
    hZ->SetLineColor(kBlue + 1);
    hZ->Draw("hist");
    c.cd(3);
    hNtrk->SetLineColor(kBlack);
    hNtrk->SetMaximum(1.15 * std::max(hNtrk->GetMaximum(), hNtpc->GetMaximum()));
    hNtrk->Draw("hist");
    hNtpc->SetLineColor(kRed);
    hNtpc->Draw("hist same");
    TLegend lg(0.5, 0.75, 0.93, 0.92);
    lg.AddEntry(hNtrk, "all tracks", "l");
    lg.AddEntry(hNtpc, "TPC-only tracks", "l");
    lg.Draw();
    hNtrk->GetXaxis()->SetTitle("tracks per event");
    save(c, "fig_event_overview");
  }

  // ================================================================ per-channel sizes and pre-cut distributions
  struct Ch
  {
    std::string tag, M, label;
  };
  std::vector<Ch> chans = {{"Ks", "K_S0", "Ks"}, {"Lambda", "Lambda0", "Lambda"}, {"D0", "D0", "Dzero"}};
  std::ofstream tsz(note + "/tables/tab_sizes.tex");
  tsz << "% generated by scripts/make_numbers.C\n";
  for (auto &c : chans)
  {
    TFile f(fn(c.tag).c_str());
    TTree *t = (TTree *) f.Get("DecayTree");
    double n = t->GetEntries(), zb = t->GetZipBytes(), tb = t->GetTotBytes();
    int nbr = t->GetListOfBranches()->GetEntries();
    num("nCand" + c.label, I(n));
    num("candPerEv" + c.label, F(n / nEv, 2));
    num("kBPerEv" + c.label, F(zb / nEv / 1024., 1));
    num("BPerCand" + c.label, I(zb / n));
    num("nBranch" + c.label, I(nbr));
    tsz << c.tag << " & " << I(n) << " & " << F(n / nEv, 2) << " & " << nbr << " & " << I(tb / n) << " & " << I(zb / n) << " & "
        << F(zb / nEv / 1024., 1) << " & " << F(zb / 1048576., 0) << " \\\\\n";

    // pre-cut variable distributions (as stored = after reconstruction pre-cuts)
    const std::string M = c.M;
    TCanvas cv(("cv" + c.tag).c_str(), "", 1500, 1000);
    cv.Divide(3, 2);
    struct V
    {
      std::string expr, xt;
      int nb;
      double lo, hi, cut;
      bool logy;
    };
    double dcaCut = c.tag == "D0" ? 0.05 : 0.5, diraCut = c.tag == "D0" ? 0.0 : 0.8;
    std::vector<V> vars = {
        {M + "_DIRA", "DIRA (cos pointing angle)", 100, -1, 1, diraCut, true},
        {"track_1_track_2_DCA", "DCA(daughter 1, daughter 2) (cm)", 100, 0, c.tag == "D0" ? 0.06 : 0.6, dcaCut, true},
        {M + "_chi2/" + M + "_nDoF", "vertex #chi^{2}/ndf", 100, 0, 25, 20, true},
        {M + "_decayLength", "decay length (cm)", 100, -5, c.tag == "D0" ? 5.0 : 40.0, -999, true},
        {"min(abs(track_1_PV_DCA_xy),abs(track_2_PV_DCA_xy))", "min daughter DCA_{xy} to PV (cm)", 100, 0, 2, -999, true},
        {M + "_pT", "p_{T} (GeV/c)", 100, 0, 5, -999, true}};
    for (size_t i = 0; i < vars.size(); ++i)
    {
      cv.cd(i + 1);
      if (vars[i].logy) gPad->SetLogy();
      std::string hn = Form("h_%s_%zu", c.tag.c_str(), i);
      t->Draw(Form("%s>>%s(%d,%f,%f)", vars[i].expr.c_str(), hn.c_str(), vars[i].nb, vars[i].lo, vars[i].hi), "", "goff");
      TH1F *h = (TH1F *) gDirectory->Get(hn.c_str());
      h->SetLineColor(kBlue + 1);
      h->GetXaxis()->SetTitle(vars[i].xt.c_str());
      h->GetYaxis()->SetTitle("candidates");
      h->Draw("hist");
      if (vars[i].cut > -998)
      {
        gPad->Update();
        TLine *ln = new TLine(vars[i].cut, std::pow(10, gPad->GetUymin()), vars[i].cut, std::pow(10, gPad->GetUymax()));
        ln->SetLineColor(kRed);
        ln->SetLineStyle(2);
        ln->SetLineWidth(2);
        ln->Draw();
      }
    }
    save(cv, "fig_precut_" + c.tag);
  }
  tsz.close();

  // ================================================================ V0 mass peaks (default, crossing-matched)
  struct Cfg
  {
    std::string tag, M, label, xt;
    double m0, s0, lo, hi, flo, fhi;
    int nb;
  };
  std::vector<Cfg> v0s = {{"Ks", "K_S0", "Ks", "m(#pi^{+}#pi^{-}) (GeV/c^{2})", 0.497611, 0.005, 0.40, 0.60, 0.43, 0.57, 100},
                          {"Lambda", "Lambda0", "Lambda", "m(p#pi) (GeV/c^{2})", 1.115683, 0.002, 1.08, 1.20, 1.09, 1.16, 120}};
  std::map<std::string, Fit> ref;
  for (auto &c : v0s)
  {
    auto mt = readCands(fn(c.tag), c.M);
    auto dup = dupFlags(mt);
    long nd = 0, ndq = 0, nq = 0, nviol = 0, nTT = 0;
    for (size_t i = 0; i < mt.size(); ++i)
    {
      if (dup[i]) { ++nd; if (qa(mt[i])) ++ndq; }
      if (qa(mt[i])) ++nq;
      if (mt[i].x1 != mt[i].x2 || mt[i].x1 != mt[i].xpv) ++nviol;
      if (mt[i].t1 == 'T' && mt[i].t2 == 'T') ++nTT;
    }
    num("dupFrac" + c.label, F(100. * nd / mt.size(), 2));
    num("dupQA" + c.label, I(ndq));
    num("nQA" + c.label, I(nq));
    num("nViolX" + c.label, I(nviol));
    num("checkX" + c.label, nviol == 0 ? "PASS" : "FAIL");
    num("nTT" + c.label, I(nTT));

    TH1F *hpre = new TH1F(("hpre" + c.tag).c_str(), "", c.nb, c.lo, c.hi);
    TH1F *hqa = new TH1F(("hqa" + c.tag).c_str(), "", c.nb, c.lo, c.hi);
    TH1F *hp = new TH1F(("hp" + c.tag).c_str(), "", c.nb, c.lo, c.hi), *ha = new TH1F(("ha" + c.tag).c_str(), "", c.nb, c.lo, c.hi);
    std::map<std::string, TH1F *> hpt;
    const std::vector<std::string> ptypes = {"II", "IM", "MM", "IT", "MT", "TT"};
    for (auto &p : ptypes) hpt[p] = new TH1F(("hpt" + c.tag + p).c_str(), "", c.nb, c.lo, c.hi);
    for (size_t i = 0; i < mt.size(); ++i)
    {
      if (dup[i]) continue;
      hpre->Fill(mt[i].mass);
      if (!qa(mt[i])) continue;
      hqa->Fill(mt[i].mass);
      hpt[pairType(mt[i].t1, mt[i].t2)]->Fill(mt[i].mass);
      if (c.tag == "Lambda") (mt[i].q > 0 ? hp : ha)->Fill(mt[i].mass);
    }
    Fit r = fitPeak(hqa, c.m0, c.s0, c.flo, c.fhi, false);
    ref[c.tag] = r;
    num("mean" + c.label, F(1e3 * r.mean, 1));
    num("emean" + c.label, F(1e3 * r.emean, 1));
    num("sigma" + c.label, F(1e3 * r.sigma, 1));
    num("esigma" + c.label, F(1e3 * r.esigma, 1));
    num("yield" + c.label, I(r.yield));
    num("eyield" + c.label, I(r.eyield));
    num("SB" + c.label, F(r.yield * 0.9973 / std::max(1.0, r.bkg3), 2));
    num("massShift" + c.label, F(1e3 * (r.mean - c.m0), 1));
    num("yieldPerKev" + c.label, F(1e3 * r.yield / nEv, 1));

    // pair-type yields with shape fixed to the inclusive fit
    std::ofstream tp(note + "/tables/tab_pairtype_" + c.tag + ".tex");
    tp << "% generated by scripts/make_numbers.C\n";
    for (auto &p : ptypes)
    {
      Fit rp = fitPeak(hpt[p], r.mean, r.sigma, c.flo, c.fhi, true, kBlue);
      tp << p << " & " << I(hpt[p]->GetEntries()) << " & " << I(rp.yield) << " $\\pm$ " << I(rp.eyield) << " & "
         << F(100 * rp.yield / r.yield, 1) << " \\\\\n";
    }
    if (c.tag == "Lambda")
    {
      Fit rp = fitPeak(hp, r.mean, r.sigma, c.flo, c.fhi, true), ra = fitPeak(ha, r.mean, r.sigma, c.flo, c.fhi, true);
      num("yieldLambdaP", I(rp.yield));
      num("eyieldLambdaP", I(rp.eyield));
      num("yieldLambdaA", I(ra.yield));
      num("eyieldLambdaA", I(ra.eyield));
    }

    TCanvas cm(("cm" + c.tag).c_str(), "", 1200, 550);
    cm.Divide(2, 1);
    cm.cd(1);
    drawMass(hpre, c.xt);
    auto l = lat();
    l.DrawLatex(0.42, 0.26, "reconstruction pre-cuts only");
    l.DrawLatex(0.42, 0.20, Form("%s candidates", I(hpre->GetEntries()).c_str()));
    cm.cd(2);
    drawMass(hqa, c.xt);
    l.DrawLatex(0.45, 0.88, "+ offline QA cuts");
    l.DrawLatex(0.45, 0.82, Form("#mu = %.1f #pm %.1f MeV", 1e3 * r.mean, 1e3 * r.emean));
    l.DrawLatex(0.45, 0.76, Form("#sigma = %.1f #pm %.1f MeV", 1e3 * r.sigma, 1e3 * r.esigma));
    l.DrawLatex(0.45, 0.70, Form("N_{sig} = %.0f #pm %.0f", r.yield, r.eyield));
    save(cm, "fig_mass_" + c.tag);

    TCanvas cp(("cp" + c.tag).c_str(), "", 1500, 1000);
    cp.Divide(3, 2);
    for (size_t i = 0; i < ptypes.size(); ++i)
    {
      cp.cd(i + 1);
      drawMass(hpt[ptypes[i]], c.xt);
      l.DrawLatex(0.5, 0.88, ("daughters " + ptypes[i]).c_str());
    }
    save(cp, "fig_pairtype_" + c.tag);
  }

  // ================================================================ crossing study: no crossing requirement
  const std::vector<std::string> classes = {"A", "B", "C", "D", "E", "F"};
  auto cls = [](const Cand &k) -> std::string {
    bool T1 = k.t1 == 'T', T2 = k.t2 == 'T';
    if (T1 && T2) return "F";
    if (T1 || T2) return ((T1 ? k.x2 : k.x1) == 0) ? "E" : "D";
    if (k.x1 == k.x2) return k.x1 == k.xpv ? "A" : "B";
    return "C";
  };
  for (auto &c : v0s)
  {
    auto nx = readCands(fn(c.tag + "_noXmatch"), c.M);
    num("nCandNoX" + c.label, I(nx.size()));
    num("candPerEvNoX" + c.label, F(nx.size() / nEv, 1));
    Fit r0 = ref[c.tag];
    std::map<std::string, TH1F *> h;
    for (auto &k : classes) h[k] = new TH1F(("hc" + c.tag + k).c_str(), "", c.nb, c.lo, c.hi);
    TH1F *hdXp = new TH1F(("hdXp" + c.tag).c_str(), "", 41, -20.5, 20.5), *hdXs = new TH1F(("hdXs" + c.tag).c_str(), "", 41, -20.5, 20.5);
    double nPk0 = 0, nSb0 = 0, nPkN = 0, nSbN = 0;  // all |dX| (no histogram range)
    // C class split by |dX|
    TH1F *hC1 = new TH1F(("hC1" + c.tag).c_str(), "", c.nb, c.lo, c.hi), *hC2 = new TH1F(("hC2" + c.tag).c_str(), "", c.nb, c.lo, c.hi);
    for (auto &k : nx)
    {
      if (!qa(k)) continue;
      std::string cl = cls(k);
      h[cl]->Fill(k.mass);
      if (cl == "C") (std::abs(k.x1 - k.x2) <= 5 ? hC1 : hC2)->Fill(k.mass);
      if (k.t1 != 'T' && k.t2 != 'T')
      {
        double d = std::fabs(k.mass - r0.mean) / r0.sigma;
        if (d < 2) { hdXp->Fill(k.x1 - k.x2); (k.x1 == k.x2 ? nPk0 : nPkN) += 1; }
        if (d > 5 && d < 7) { hdXs->Fill(k.x1 - k.x2); (k.x1 == k.x2 ? nSb0 : nSbN) += 1; }  // same total width (2x2 sigma) as peak window
      }
    }
    std::ofstream tc(note + "/tables/tab_xclass_" + c.tag + ".tex");
    tc << "% generated by scripts/make_numbers.C\n";
    std::map<std::string, std::string> desc = {{"A", "$X_1=X_2=X_{\\rm PV}$ (both Si)"},
                                               {"B", "$X_1=X_2\\neq X_{\\rm PV}$ (both Si)"},
                                               {"C", "$X_1\\neq X_2$ (both Si)"},
                                               {"D", "TPC-only + Si with $X\\neq0$"},
                                               {"E", "TPC-only + Si with $X=0$"},
                                               {"F", "both TPC-only"}};
    for (auto &k : classes)
    {
      Fit rk = fitPeak(h[k], r0.mean, r0.sigma, c.flo, c.fhi, true, kBlue);
      tc << k << " & " << desc[k] << " & " << I(h[k]->GetEntries()) << " & " << I(rk.yield) << " $\\pm$ " << I(rk.eyield) << " & "
         << F(rk.eyield > 0 ? rk.yield / rk.eyield : 0, 1) << " \\\\\n";
      num("xcl" + k + c.label, I(rk.yield));
      num("excl" + k + c.label, I(rk.eyield));
    }
    Fit rc1 = fitPeak(hC1, r0.mean, r0.sigma, c.flo, c.fhi, true, kBlue), rc2 = fitPeak(hC2, r0.mean, r0.sigma, c.flo, c.fhi, true, kBlue);
    num("xclCsmall" + c.label, I(rc1.yield));
    num("exclCsmall" + c.label, I(rc1.eyield));
    num("xclClarge" + c.label, I(rc2.yield));
    num("exclClarge" + c.label, I(rc2.eyield));

    TCanvas cc(("cc" + c.tag).c_str(), "", 1500, 1000);
    cc.Divide(3, 2);
    auto l = lat();
    for (size_t i = 0; i < classes.size(); ++i)
    {
      cc.cd(i + 1);
      drawMass(h[classes[i]], c.xt);
      l.DrawLatex(0.45, 0.88, ("class " + classes[i]).c_str());
    }
    save(cc, "fig_xclass_" + c.tag);

    TCanvas cd(("cd" + c.tag).c_str(), "", 700, 550);
    TH1F *hsig = (TH1F *) hdXp->Clone(("hsig" + c.tag).c_str());
    hsig->Add(hdXs, -1);
    hdXp->SetLineColor(kBlack);
    hdXs->SetLineColor(kRed);
    hdXs->SetLineWidth(2);
    hdXp->SetLineWidth(2);
    hsig->SetMarkerStyle(20);
    hsig->SetMarkerColor(kBlue + 1);
    hsig->SetLineColor(kBlue + 1);
    double mx = std::max(hdXp->GetMaximum(), hsig->GetMaximum());
    hdXp->SetMaximum(1.3 * mx);
    hdXp->SetMinimum(std::min(0.0, 1.2 * hsig->GetMinimum()));
    hdXp->GetXaxis()->SetTitle("X_{1} - X_{2} (bunch crossings)");
    hdXp->GetYaxis()->SetTitle("candidates");
    hdXp->Draw("hist");
    hdXs->Draw("hist same");
    hsig->Draw("E same");
    TLegend lg(0.55, 0.72, 0.95, 0.93);
    lg.AddEntry(hdXp, "peak window (|m-#mu|<2#sigma)", "l");
    lg.AddEntry(hdXs, "sidebands (5-7#sigma)", "l");
    lg.AddEntry(hsig, "peak - sidebands", "lep");
    lg.Draw();
    save(cd, "fig_dX_" + c.tag);
    double s0 = hsig->GetBinContent(hsig->FindBin(0)), e0 = hsig->GetBinError(hsig->FindBin(0));
    double sall = 0, eall2 = 0;
    for (int b = 1; b <= hsig->GetNbinsX(); ++b) { sall += hsig->GetBinContent(b); eall2 += std::pow(hsig->GetBinError(b), 2); }
    num("dXzeroAll" + c.label, F(nPk0 - nSb0, 0));
    num("edXzeroAll" + c.label, F(std::sqrt(nPk0 + nSb0), 0));
    num("dXnonzeroAll" + c.label, F(nPkN - nSbN, 0));
    num("edXnonzeroAll" + c.label, F(std::sqrt(nPkN + nSbN), 0));
    num("dXnonzeroFrac" + c.label, F(100 * (nPkN - nSbN) / std::max(1.0, nPk0 - nSb0 + nPkN - nSbN), 0));
    num("dXzero" + c.label, F(s0, 0));
    num("edXzero" + c.label, F(e0, 0));
    num("dXnonzero" + c.label, F(sall - s0, 0));
    num("edXnonzero" + c.label, F(std::sqrt(std::max(0.0, eall2 - e0 * e0)), 0));
  }

  // ================================================================ D0
  {
    TFile f(fn("D0").c_str());
    TTree *t = (TTree *) f.Get("DecayTree");
    TTreeReader r(t);
    TTreeReaderValue<float> mass(r, "D0_mass"), dira(r, "D0_DIRA"), dl(r, "D0_decayLength"), dle(r, "D0_decayLengthErr"),
        pt(r, "D0_pT"), dca(r, "track_1_track_2_DCA"), p1x(r, "track_1_px"), p1y(r, "track_1_py"), p1z(r, "track_1_pz"),
        p2x(r, "track_2_px"), p2y(r, "track_2_py"), p2z(r, "track_2_pz"), d1(r, "track_1_PV_DCA_xy"), d2(r, "track_2_PV_DCA_xy");
    TTreeReaderValue<int> pdg1(r, "track_1_PDG_ID");
    TTreeReaderValue<unsigned> i1(r, "track_1_INTT_nHits"), i2(r, "track_2_INTT_nHits");
    TH1F *hm = new TH1F("hD0m", "", 100, 1.6, 2.1), *hs = new TH1F("hD0s", "", 100, 1.6, 2.1), *hsel = new TH1F("hD0sel", "", 50, 1.6, 2.1);
    const double mK = 0.493677, mpi = 0.13957;
    long nsel = 0;
    while (r.Next())
    {
      hm->Fill(*mass);
      // swapped hypothesis: the track that was the kaon becomes the pion and vice versa
      bool t1K = std::abs(*pdg1) == 321;
      double ma = t1K ? mpi : mK, mb = t1K ? mK : mpi;
      double pa2 = (*p1x) * (*p1x) + (*p1y) * (*p1y) + (*p1z) * (*p1z), pb2 = (*p2x) * (*p2x) + (*p2y) * (*p2y) + (*p2z) * (*p2z);
      double E = std::sqrt(pa2 + ma * ma) + std::sqrt(pb2 + mb * mb);
      double px = *p1x + *p2x, py = *p1y + *p2y, pz = *p1z + *p2z;
      hs->Fill(std::sqrt(std::max(0.0, E * E - px * px - py * py - pz * pz)));
      bool sel = *dira > 0.95 && *dle > 0 && *dl / *dle > 3 && *dca < 0.01 && std::fabs(*d1) > 0.005 && std::fabs(*d2) > 0.005 && *i1 > 0 && *i2 > 0;
      if (sel) { hsel->Fill(*mass); ++nsel; }
    }
    num("nDzeroSel", I(nsel));
    TCanvas c("cD0", "", 1200, 550);
    c.Divide(2, 1);
    c.cd(1);
    hm->SetLineColor(kBlack);
    hm->GetXaxis()->SetTitle("m(K#pi) (GeV/c^{2})");
    hm->GetYaxis()->SetTitle("candidates / 5 MeV");
    hm->SetMinimum(0);
    hm->Draw("hist");
    hs->SetLineColor(kRed);
    hs->Draw("hist same");
    TLegend lg(0.5, 0.78, 0.95, 0.93);
    lg.AddEntry(hm, "stored hypothesis", "l");
    lg.AddEntry(hs, "swapped K#leftrightarrow#pi", "l");
    lg.Draw();
    c.cd(2);
    drawMass(hsel, "m(K#pi) (GeV/c^{2})");
    auto l = lat();
    l.DrawLatex(0.2, 0.25, "DIRA>0.95, L/#sigma_{L}>3, DCA_{12}<100 #mum");
    l.DrawLatex(0.2, 0.19, "DCA_{xy}^{trk}>50 #mum, both INTT");
    save(c, "fig_mass_D0");
  }

  fclose(gLog);
  gMac.close();
  printf("make_numbers done: %s\n", (note + "/logs/numbers.txt").c_str());
}
