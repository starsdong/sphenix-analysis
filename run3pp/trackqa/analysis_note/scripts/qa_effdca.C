// qa_effdca.C : silicon-matching efficiency and DCA vs p, eta, phi and vs crossing / vertex z / multiplicity.
//   root -l -b -q 'scripts/qa_effdca.C("<note root>")'
// Input : <note>/../ntuples/TrackQA_run79509_seg000-099.root (100 segments, hadd)
// Output: figures/fig_eff_*.{pdf,png}, figures/fig_dcav_*.{pdf,png}, tables/effdca_numbers.tex,
//         tables/effdca_defaults.tex, tables/tab_eff_pt.tex, logs/effdca_numbers.txt
//
// Silicon efficiency  = N(ntpc>=20 && nintt>=1 && nmvtx>=2) / N(ntpc>=20), all tracks of SvtxTrackMap.
// DCA sample          = ntpc>=20 && nintt>=1 && nmvtx>=2 && vtxMatch==1 (track used in its vertex fit).
// eta, phi, crossing, vz, multiplicity dependences: pT > 1 GeV/c.
// Gaussian core mean/sigma per x bin: iterative fit in mean +- 2 sigma on 10-um bins (same as dca_study.C).
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TLine.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TEfficiency.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TMath.h>
#include <cmath>
#include <fstream>
#include <iostream>
#include <regex>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

namespace ed
{
  std::string note;
  std::ofstream tex, txt;
  int nFail = 0, nCheck = 0;
  std::string tagText;
  void num(const std::string &n, const std::string &v)
  {
    tex << "\\renewcommand{\\" << n << "}{" << v << "}\n";
    txt << n << " = " << v << "\n";
  }
  void numf(const std::string &n, double v, int p)
  {
    std::ostringstream s;
    s.setf(std::ios::fixed);
    s.precision(p);
    if (std::fabs(v) < 0.5 * std::pow(10., -p)) v = 0;  // no "-0"
    s << v;
    num(n, s.str());
  }
  void check(const std::string &n, bool ok, const std::string &what)
  {
    num(n, ok ? "PASS" : "FAIL");
    txt << "  CHECK " << n << ": " << what << " -> " << (ok ? "PASS" : "FAIL") << "\n";
    ++nCheck;
    if (!ok) ++nFail;
  }
  std::string fmtInt(long v)
  {
    std::string s = std::to_string(v), o;
    int n = 0;
    for (int i = s.size() - 1; i >= 0; --i)
    {
      o.insert(o.begin(), s[i]);
      if (++n % 3 == 0 && i > 0) o.insert(o.begin(), ',');
    }
    return o;
  }
  void save(TCanvas *c, const std::string &n)
  {
    c->SaveAs((note + "/figures/" + n + ".pdf").c_str());
    c->SaveAs((note + "/figures/" + n + ".png").c_str());
  }
  void ndcLabel(const std::string &t, double x = 0.2, double y = 0.935, double sz = 0.05)
  {
    TLatex l;
    l.SetNDC();
    l.SetTextFont(42);
    l.SetTextSize(sz);
    l.DrawLatex(x, y, t.c_str());
  }
  void pad(bool right = false)
  {
    gPad->SetTopMargin(0.1);
    gPad->SetRightMargin(right ? 0.16 : 0.04);
    gPad->SetLeftMargin(0.2);
    gPad->SetBottomMargin(0.16);
  }
  struct Fit
  {
    double mean = 0, emean = 0, sig = -1, esig = 0;
  };
  Fit core(TH1 *h)
  {
    Fit r;
    if (h->GetEntries() < 100) return r;
    double m = h->GetBinCenter(h->GetMaximumBin());
    double s = std::min(h->GetRMS(), 0.05);
    TF1 f("g", "gaus", m - 2 * s, m + 2 * s);
    for (int it = 0; it < 4; ++it)
    {
      f.SetRange(m - 2 * s, m + 2 * s);
      f.SetParameters(h->GetMaximum(), m, s);
      h->Fit(&f, "QNR0");
      m = f.GetParameter(1);
      s = std::fabs(f.GetParameter(2));
    }
    r.mean = m;
    r.emean = f.GetParError(1);
    r.sig = s;
    r.esig = f.GetParError(2);
    return r;
  }
  // one DCA-vs-x study: 2D histograms (10-um y bins) for xy and z
  struct DV
  {
    std::string name, xtitle, sel;
    bool logx;
    TH2D *h[2];
    TGraphErrors *gS[2], *gM[2];
  };
  DV makeDV(const std::string &name, const std::string &xt, const std::string &sel, int nb, const double *edges, bool logx)
  {
    DV d;
    d.name = name;
    d.xtitle = xt;
    d.sel = sel;
    d.logx = logx;
    for (int a = 0; a < 2; ++a)
    {
      d.h[a] = new TH2D(Form("h2_%s_%d", name.c_str(), a), "", nb, edges, 1000, -0.5, 0.5);
      d.gS[a] = new TGraphErrors();
      d.gM[a] = new TGraphErrors();
    }
    return d;
  }
  void fitDV(DV &d)
  {
    for (int a = 0; a < 2; ++a)
      for (int b = 1; b <= d.h[a]->GetNbinsX(); ++b)
      {
        TH1D *p = d.h[a]->ProjectionY(Form("py_%s_%d_%d", d.name.c_str(), a, b), b, b);
        Fit r = core(p);
        if (r.sig < 0) continue;
        double xc = d.h[a]->GetXaxis()->GetBinCenter(b), ex = 0.5 * d.h[a]->GetXaxis()->GetBinWidth(b);
        int n = d.gS[a]->GetN();
        d.gS[a]->SetPoint(n, xc, r.sig * 1e4);
        d.gS[a]->SetPointError(n, ex, r.esig * 1e4);
        d.gM[a]->SetPoint(n, xc, r.mean * 1e4);
        d.gM[a]->SetPointError(n, ex, r.emean * 1e4);
      }
  }
  void minmax(TGraph *g, double &mn, double &mx)
  {
    mn = 1e9;
    mx = -1e9;
    for (int i = 0; i < g->GetN(); ++i)
    {
      mn = std::min(mn, g->GetY()[i]);
      mx = std::max(mx, g->GetY()[i]);
    }
  }
  // 4 panels: 2D xy, 2D z (log z, +-0.1 cm), core sigma, core mean
  void drawDV(DV &d, double sMax, double mRange)
  {
    auto *c = new TCanvas(("c" + d.name).c_str(), "", 1800, 470);
    c->Divide(4, 1, 0.002, 0.002);
    for (int a = 0; a < 2; ++a)
    {
      c->cd(a + 1);
      pad(true);
      if (d.logx) gPad->SetLogx();
      gPad->SetLogz();
      TH2D *hc = (TH2D *) d.h[a]->Clone(Form("%s_draw", d.h[a]->GetName()));
      hc->RebinY(4);
      hc->GetYaxis()->SetRangeUser(-0.1, 0.0999);
      hc->GetXaxis()->SetTitle(d.xtitle.c_str());
      hc->GetYaxis()->SetTitle(a ? "DCA_{z} (cm)" : "DCA_{xy} (cm)");
      hc->GetXaxis()->SetMoreLogLabels();
      hc->Draw("colz");
      ndcLabel(a ? "DCA_{z}" : (tagText + ", " + d.sel), 0.2, 0.935, a ? 0.05 : 0.042);
    }
    const char *yt[2] = {"core #sigma (#mum)", "core mean (#mum)"};
    for (int q = 0; q < 2; ++q)
    {
      c->cd(q + 3);
      pad();
      if (d.logx) gPad->SetLogx();
      double x0 = d.h[0]->GetXaxis()->GetXmin(), x1 = d.h[0]->GetXaxis()->GetXmax();
      TH1F *fr = gPad->DrawFrame(x0, q ? -mRange : 0, x1, q ? mRange : sMax);
      fr->SetXTitle(d.xtitle.c_str());
      fr->SetYTitle(yt[q]);
      fr->GetXaxis()->SetMoreLogLabels();
      if (q)
      {
        TLine ln;
        ln.SetLineStyle(2);
        ln.SetLineColor(kGray + 1);
        ln.DrawLine(x0, 0, x1, 0);
      }
      auto *lg = new TLegend(0.25, 0.78, 0.95, 0.89);
      lg->SetBorderSize(0);
      lg->SetTextSize(0.05);
      lg->SetNColumns(2);
      for (int a = 0; a < 2; ++a)
      {
        TGraphErrors *g = q ? d.gM[a] : d.gS[a];
        g->SetMarkerStyle(a ? 21 : 20);
        g->SetMarkerColor(a ? kRed + 1 : kBlack);
        g->SetLineColor(a ? kRed + 1 : kBlack);
        g->SetMarkerSize(1.0);
        g->Draw("p same");
        lg->AddEntry(g, a ? "DCA_{z}" : "DCA_{xy}", "p");
      }
      lg->Draw();
    }
    save(c, "fig_dcav_" + d.name);
  }
  // efficiency vs x
  struct EV
  {
    std::string name, xtitle;
    TH1D *pass, *tot;
    TGraphAsymmErrors *g = nullptr;
  };
  EV makeEV(const std::string &name, const std::string &xt, int nb, const double *edges)
  {
    EV e;
    e.name = name;
    e.xtitle = xt;
    e.pass = new TH1D(("ep_" + name).c_str(), "", nb, edges);
    e.tot = new TH1D(("et_" + name).c_str(), "", nb, edges);
    return e;
  }
  EV makeEVu(const std::string &name, const std::string &xt, int nb, double lo, double hi)
  {
    std::vector<double> ed(nb + 1);
    for (int i = 0; i <= nb; ++i) ed[i] = lo + (hi - lo) * i / nb;
    return makeEV(name, xt, nb, ed.data());
  }
  void effGraph(EV &e)
  {
    TEfficiency te(*e.pass, *e.tot);
    e.g = new TGraphAsymmErrors();
    for (int b = 1; b <= e.tot->GetNbinsX(); ++b)
    {
      if (e.tot->GetBinContent(b) < 20) continue;
      int n = e.g->GetN();
      double xc = e.tot->GetBinCenter(b), ex = 0.5 * e.tot->GetBinWidth(b);
      e.g->SetPoint(n, xc, te.GetEfficiency(b));
      e.g->SetPointError(n, ex, ex, te.GetEfficiencyErrorLow(b), te.GetEfficiencyErrorUp(b));
    }
  }
  void drawEff(EV &e, bool logx, double ylo = 0.0, double yhi = 1.05, const std::string &txt = "")
  {
    gPad->SetLogx(logx);
    TH1F *fr = gPad->DrawFrame(e.tot->GetXaxis()->GetXmin(), ylo, e.tot->GetXaxis()->GetXmax(), yhi);
    fr->SetXTitle(e.xtitle.c_str());
    fr->SetYTitle("silicon efficiency");
    fr->GetXaxis()->SetMoreLogLabels();
    e.g->SetMarkerStyle(20);
    e.g->SetMarkerSize(0.9);
    e.g->Draw("p same");
    if (!txt.empty()) ndcLabel(txt, 0.26, 0.24, 0.05);
  }
}  // namespace ed
using namespace ed;

void qa_effdca(const char *noteRoot = ".")
{
  note = noteRoot;
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetPalette(kBird);
  TH1::SetDefaultSumw2();
  tex.open(note + "/tables/effdca_numbers.tex");
  txt.open(note + "/logs/effdca_numbers.txt");
  tex << "% generated by scripts/qa_effdca.C -- do not edit\n";

  TFile *f = TFile::Open((note + "/../ntuples/TrackQA_run79509_seg000-099.root").c_str());
  TTree *te = (TTree *) f->Get("event");
  TTree *tt = (TTree *) f->Get("track");
  const long nEv = te->GetEntries(), nTr = tt->GetEntries();
  tagText = "run3pp 79509 seg. 0-99";

  // ---- event map: event number -> entry
  Int_t eEvt, eNTrk, eBunch, eNVtx0, eNSilX0;
  Float_t eVz;
  te->SetBranchAddress("event", &eEvt);
  te->SetBranchAddress("nTrk", &eNTrk);
  te->SetBranchAddress("GL1_bunch", &eBunch);
  te->SetBranchAddress("nVtxX0", &eNVtx0);
  te->SetBranchAddress("nTrkSilX0", &eNSilX0);
  te->SetBranchAddress("vz", &eVz);
  std::unordered_map<int, long> evIdx;
  std::vector<int> vNTrk(nEv), vBunch(nEv), vNVtx0(nEv);
  std::vector<float> vVz(nEv);
  long nDup = 0;
  int bunchMin = 1e9, bunchMax = -1;
  for (long i = 0; i < nEv; ++i)
  {
    te->GetEntry(i);
    if (evIdx.count(eEvt)) ++nDup;
    else evIdx[eEvt] = i;
    vNTrk[i] = eNSilX0;  // multiplicity of the triggered crossing (silicon tracks at crossing 0)
    vBunch[i] = eBunch;
    vNVtx0[i] = eNVtx0;
    vVz[i] = eVz;
    bunchMin = std::min(bunchMin, eBunch);
    bunchMax = std::max(bunchMax, eBunch);
  }
  check("chkEDup", nDup == 0, "event numbers unique in the merged 100-segment event tree");
  num("eBunchMin", std::to_string(bunchMin));
  num("eBunchMax", std::to_string(bunchMax));

  // ---- binning
  const double ptE[] = {0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 1.0, 1.2, 1.5, 2.0, 2.5, 3.0, 4.0, 6.0, 10.0};
  const int nPtE = sizeof(ptE) / sizeof(double) - 1;
  const double pE[] = {0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 1.0, 1.2, 1.5, 2.0, 2.5, 3.0, 4.0, 6.0, 10.0};
  const int nPE = sizeof(pE) / sizeof(double) - 1;
  std::vector<double> etaE, phiE, xE, vzE, ntE, entE, buE;
  for (int i = 0; i <= 24; ++i) etaE.push_back(-1.2 + 0.1 * i);  // DCA
  std::vector<double> etaEff;
  for (int i = 0; i <= 30; ++i) etaEff.push_back(-1.5 + 0.1 * i);  // efficiency acceptance |eta|<1.5
  for (int i = 0; i <= 72; ++i) phiE.push_back(-TMath::Pi() + i * TMath::Pi() / 36);
  for (int i = 0; i <= 32; ++i) xE.push_back(-120 + 20 * i);  // track crossing -120..520
  for (int i = 0; i <= 20; ++i) vzE.push_back(-25 + 2.5 * i);
  for (double v : {0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 7.5, 9.5, 12.5, 16.5, 24.5, 39.5, 80.5}) ntE.push_back(v);
  for (int i = 0; i <= 50; ++i) entE.push_back(i * 2.);  // event entry / 1000 (time order)
  for (int i = 0; i <= 120; ++i) buE.push_back(i - 0.5);

  // efficiency
  EV eP = makeEV("pt", "#it{p}_{T} (GeV/#it{c})", nPtE, ptE);
  EV eEta = makeEV("eta", "#eta", etaEff.size() - 1, etaEff.data());
  EV ePhi = makeEV("phi", "#phi (rad)", phiE.size() - 1, phiE.data());
  EV eBu = makeEV("bunch", "GL1 bunch number", buE.size() - 1, buE.data());
  EV eVzE = makeEV("vz", "event #it{v}_{z} (crossing-0 vertex, cm)", vzE.size() - 1, vzE.data());
  EV eNt = makeEV("ntrk", "#it{N}_{Si} at crossing 0 (event)", ntE.size() - 1, ntE.data());
  EV eEnt = makeEV("entry", "event entry / 1000 (#approx segment)", entE.size() - 1, entE.data());
  EV ePz = makeEV("pcaz", "track #it{z} at PCA (cm)", vzE.size() - 1, vzE.data());
  TH2D *h2pass = new TH2D("h2pass", ";#eta;#phi (rad)", 30, -1.5, 1.5, 72, -TMath::Pi(), TMath::Pi());
  TH2D *h2tot = (TH2D *) h2pass->Clone("h2tot");
  TH2D *h2ptPass = new TH2D("h2ptPass", ";#eta;#it{p}_{T} (GeV/#it{c})", 30, -1.5, 1.5, nPtE, ptE);
  TH2D *h2ptTot = (TH2D *) h2ptPass->Clone("h2ptTot");

  // DCA
  DV dP = makeDV("p", "#it{p} (GeV/#it{c})", "all #it{p}_{T}", nPE, pE, true);
  DV dEta = makeDV("eta", "#eta", "#it{p}_{T} > 1 GeV/#it{c}", etaE.size() - 1, etaE.data(), false);
  DV dPhi = makeDV("phi", "#phi (rad)", "#it{p}_{T} > 1 GeV/#it{c}", phiE.size() - 1, phiE.data(), false);
  DV dX = makeDV("crossing", "track beam crossing", "#it{p}_{T} > 1 GeV/#it{c}", xE.size() - 1, xE.data(), false);
  DV dVz = makeDV("vz", "vertex #it{z} (cm)", "#it{p}_{T} > 1 GeV/#it{c}", vzE.size() - 1, vzE.data(), false);
  DV dNt = makeDV("ntrk", "#it{N}_{Si} in same crossing", "#it{p}_{T} > 1 GeV/#it{c}", ntE.size() - 1, ntE.data(), true);
  DV dEnt = makeDV("entry", "event entry / 1000 (#approx segment)", "#it{p}_{T} > 1 GeV/#it{c}", entE.size() - 1, entE.data(), false);

  // ---- track loop
  Int_t tEvt;
  Float_t pt, eta, phi, dcaxy, dcaz, tvz, pcaz;
  Short_t vtxMatch, crossing;
  UChar_t nmvtx, nintt, ntpc;
  tt->SetBranchStatus("*", 0);
  for (const char *b : {"event", "pt", "eta", "phi", "dcaxy", "dcaz", "tvz", "pcaz", "vtxMatch", "crossing", "nmvtx", "nintt", "ntpc"})
    tt->SetBranchStatus(b, 1);
  tt->SetBranchAddress("event", &tEvt);
  tt->SetBranchAddress("pt", &pt);
  tt->SetBranchAddress("eta", &eta);
  tt->SetBranchAddress("phi", &phi);
  tt->SetBranchAddress("dcaxy", &dcaxy);
  tt->SetBranchAddress("dcaz", &dcaz);
  tt->SetBranchAddress("tvz", &tvz);
  tt->SetBranchAddress("pcaz", &pcaz);
  tt->SetBranchAddress("vtxMatch", &vtxMatch);
  tt->SetBranchAddress("crossing", &crossing);
  tt->SetBranchAddress("nmvtx", &nmvtx);
  tt->SetBranchAddress("nintt", &nintt);
  tt->SetBranchAddress("ntpc", &ntpc);
  // ---- first pass: number of silicon tracks (MVTX or INTT cluster, any pT) per (event, crossing)
  std::unordered_map<long long, int> nSil;
  auto xkey = [](int ev, int x) { return (long long) ev * 4096LL + (x + 2048); };
  long nSilTot = 0;
  for (long i = 0; i < nTr; ++i)
  {
    tt->GetEntry(i);
    if (nmvtx + nintt > 0)
    {
      ++nSil[xkey(tEvt, crossing)];
      ++nSilTot;
    }
  }
  {
    TH1D hm("hm", "", 200, 0.5, 200.5);
    for (auto &kv : nSil) hm.Fill(kv.second);
    numf("eSilPerXing", hm.GetMean(), 1);
    numf("eXingPerEv", double(nSil.size()) / nEv, 1);
    double q = 0.5, med;
    hm.GetQuantiles(1, &med, &q);
    numf("eSilPerXingMedian", med, 0);
  }
  long nNoEv = 0, nDen = 0, nNum = 0, nDen1 = 0, nNum1 = 0, nDca = 0, nDca1 = 0, nNumMatch1 = 0;
  for (long i = 0; i < nTr; ++i)
  {
    tt->GetEntry(i);
    if (ntpc < 20) continue;
    auto it = evIdx.find(tEvt);
    if (it == evIdx.end())
    {
      ++nNoEv;
      continue;
    }
    const long ie = it->second;
    const bool pass = nintt >= 1 && nmvtx >= 2;
    const bool hi = pt > 1.0;
    // matching-efficiency acceptance (user choice): track z at PCA |z|<10 cm, |eta|<1.5; vs pT: |eta|<1.0.
    // The z-dependence plots (event vz, track z) apply only |eta|<1.5.
    const bool inEta = std::fabs(eta) < 1.5, inZ = std::fabs(pcaz) < 10;
    const bool acc = inEta && inZ;
    auto fillE = [&](EV &e, double x) {
      e.tot->Fill(x);
      if (pass) e.pass->Fill(x);
    };
    if (acc && std::fabs(eta) < 1.0)
    {
      ++nDen;
      if (pass) ++nNum;
      fillE(eP, pt);
    }
    if (acc)
    {
      h2ptTot->Fill(eta, pt);
      if (pass) h2ptPass->Fill(eta, pt);
    }
    if (hi && inEta)
    {
      if (vNVtx0[ie] > 0) fillE(eVzE, vVz[ie]);
      fillE(ePz, pcaz);
    }
    if (hi && acc)
    {
      ++nDen1;
      if (pass) ++nNum1;
      fillE(eEta, eta);
      fillE(ePhi, phi);
      fillE(eBu, vBunch[ie]);
      fillE(eNt, vNTrk[ie]);
      fillE(eEnt, ie / 1000.);
      h2tot->Fill(eta, phi);
      if (pass) h2pass->Fill(eta, phi);
    }
    // DCA sample
    if (!pass || vtxMatch != 1) continue;
    ++nDca;
    const double p = pt * std::cosh(eta);
    const double d[2] = {dcaxy, dcaz};
    for (int a = 0; a < 2; ++a) dP.h[a]->Fill(p, d[a]);
    if (!hi) continue;
    ++nDca1;
    const int nSameX = nSil[xkey(tEvt, crossing)];
    if (crossing != 0 || true) ++nNumMatch1;
    for (int a = 0; a < 2; ++a)
    {
      dEta.h[a]->Fill(eta, d[a]);
      dPhi.h[a]->Fill(phi, d[a]);
      dX.h[a]->Fill(crossing, d[a]);
      dVz.h[a]->Fill(tvz, d[a]);
      dNt.h[a]->Fill(nSameX, d[a]);
      dEnt.h[a]->Fill(ie / 1000., d[a]);
    }
  }
  check("chkENoEv", nNoEv == 0, "every track found its event in the event tree (by event number)");
  check("chkECount", std::abs(eP.tot->GetEntries() - nDen) < 0.5 && std::abs(eEta.tot->GetEntries() - nDen1) < 0.5,
        "efficiency denominators equal the counted tracks (vs pT: |z|<10,|eta|<1; pT>1: |z|<10,|eta|<1.5)");
  num("eNEvents", fmtInt(nEv));
  num("eNDen", fmtInt(nDen));
  num("eNDenHi", fmtInt(nDen1));
  num("eNDca", fmtInt(nDca));
  num("eNDcaHi", fmtInt(nDca1));
  numf("eEffAll", 100. * nNum / nDen, 1);
  numf("eEffHi", 100. * nNum1 / nDen1, 1);

  // ---- efficiency graphs + numbers
  for (EV *e : {&eP, &eEta, &ePhi, &eBu, &eVzE, &eNt, &eEnt, &ePz}) effGraph(*e);
  auto effRange = [&](EV &e, const std::string &tag) {
    double mn, mx;
    minmax(e.g, mn, mx);
    numf("eEff" + tag + "Min", 100 * mn, 1);
    numf("eEff" + tag + "Max", 100 * mx, 1);
  };
  effRange(eEta, "Eta");
  effRange(ePhi, "Phi");
  effRange(eBu, "Bunch");
  effRange(eVzE, "Vz");
  effRange(eNt, "Ntrk");
  effRange(eEnt, "Entry");
  effRange(ePz, "Pcaz");
  {
    // phi bins (pT>1) with efficiency below 70%
    std::string lst;
    int nlow = 0;
    for (int i = 0; i < ePhi.g->GetN(); ++i)
      if (ePhi.g->GetY()[i] < 0.70)
      {
        lst += (nlow ? ", " : "") + std::string(Form("%.2f", ePhi.g->GetX()[i]));
        ++nlow;
      }
    num("eLowPhiList", lst);
    num("eLowPhiN", std::to_string(nlow));
    // pcaz plateau: |z|<10
    double s1 = 0, n1 = 0;
    for (int i = 0; i < ePz.g->GetN(); ++i)
      if (std::fabs(ePz.g->GetX()[i]) < 10)
      {
        s1 += ePz.g->GetY()[i];
        ++n1;
      }
    numf("eEffPcazPlateau", 100 * s1 / n1, 1);
  }
  numf("eEffPtLow", 100 * eP.g->GetY()[0], 1);
  numf("eEffPtHigh", 100 * eP.g->GetY()[eP.g->GetN() - 1], 1);
  {
    double mn, mx;
    minmax(eP.g, mn, mx);
    numf("eEffPtMax", 100 * mx, 1);
  }
  // pT table
  {
    std::ofstream t(note + "/tables/tab_eff_pt.tex");
    t << "% generated by scripts/qa_effdca.C\n\\begin{tabular}{lrrr}\n\\toprule\n$p_T$ (GeV/$c$) & denominator & numerator & efficiency (\\%)\\\\\n\\midrule\n";
    TEfficiency tef(*eP.pass, *eP.tot);
    for (int b = 1; b <= eP.tot->GetNbinsX(); ++b)
      t << Form("%.1f--%.1f & %s & %s & %.1f $^{+%.1f}_{-%.1f}$\\\\\n", ptE[b - 1], ptE[b], fmtInt((long) eP.tot->GetBinContent(b)).c_str(),
                fmtInt((long) eP.pass->GetBinContent(b)).c_str(), 100 * tef.GetEfficiency(b), 100 * tef.GetEfficiencyErrorUp(b),
                100 * tef.GetEfficiencyErrorLow(b));
    t << "\\bottomrule\n\\end{tabular}\n";
  }

  // ---- efficiency figures
  {
    auto *c = new TCanvas("ceff1", "", 1650, 500);
    c->Divide(3, 1, 0.002, 0.002);
    c->cd(1);
    pad();
    drawEff(eP, true, 0, 1.05, "|#it{z}|<10 cm, |#eta|<1");
    ndcLabel(tagText + ", " + fmtInt(nEv) + " events", 0.2, 0.935, 0.045);
    c->cd(2);
    pad();
    drawEff(eEta, false, 0, 1.05, "#it{p}_{T}>1, |#it{z}|<10 cm");
    c->cd(3);
    pad();
    drawEff(ePhi, false, 0, 1.05, "#it{p}_{T}>1, |#it{z}|<10, |#eta|<1.5");
    save(c, "fig_eff_kin");
  }
  {
    auto *c = new TCanvas("ceff2", "", 1200, 500);
    c->Divide(2, 1, 0.002, 0.002);
    TH2D *r1 = (TH2D *) h2pass->Clone("r1");
    r1->Divide(h2pass, h2tot, 1, 1, "B");
    TH2D *r2 = (TH2D *) h2ptPass->Clone("r2");
    r2->Divide(h2ptPass, h2ptTot, 1, 1, "B");
    TH2D *rr[2] = {r1, r2};
    for (int i = 0; i < 2; ++i)
    {
      c->cd(i + 1);
      pad(true);
      if (i == 1) gPad->SetLogy();
      rr[i]->SetMinimum(0);
      rr[i]->SetMaximum(1);
      rr[i]->GetZaxis()->SetTitle("silicon efficiency");
      rr[i]->GetZaxis()->SetTitleOffset(1.1);
      rr[i]->GetYaxis()->SetMoreLogLabels();
      rr[i]->GetYaxis()->SetTitleSize(0.055);
      rr[i]->GetYaxis()->SetTitleOffset(1.3);
      rr[i]->GetXaxis()->SetTitleSize(0.055);
      rr[i]->Draw("colz");
      ndcLabel(i ? "|#it{z}|<10 cm" : tagText + ", #it{p}_{T}>1, |#it{z}|<10 cm", 0.2, 0.935, 0.045);
    }
    save(c, "fig_eff_2d");
  }
  {
    auto *c = new TCanvas("ceff3", "", 1800, 900);
    c->Divide(3, 2, 0.002, 0.002);
    EV *es[6] = {&eBu, &eVzE, &ePz, &eNt, &eEnt, nullptr};
    for (int i = 0; i < 5; ++i)
    {
      c->cd(i + 1);
      pad();
      drawEff(*es[i], false, 0, 1.05, (i == 1 || i == 2) ? "#it{p}_{T}>1, |#eta|<1.5" : "#it{p}_{T}>1, |#it{z}|<10, |#eta|<1.5");
      if (i == 0) ndcLabel(tagText, 0.2, 0.935, 0.05);
    }
    c->cd(6);
    pad();
    // track crossing of the numerator tracks (denominator crossing undefined for TPC-only tracks)
    TH1D *hx = new TH1D("hxnum", ";track beam crossing;silicon tracks (#it{p}_{T} > 1)", 64, -120, 520);
    tt->Draw("crossing>>hxnum", "ntpc>=20&&nintt>=1&&nmvtx>=2&&pt>1", "goff");
    hx->SetLineColor(kBlue + 1);
    hx->SetLineWidth(2);
    hx->Draw("hist");
    gPad->SetLogy();
    ndcLabel("numerator only", 0.5, 0.82, 0.05);
    numf("eFracNumXzero", 100. * hx->GetBinContent(hx->FindBin(0)) / hx->GetEntries(), 1);
    save(c, "fig_eff_event");
  }

  // ---- DCA studies
  for (DV *d : {&dP, &dEta, &dPhi, &dX, &dVz, &dNt, &dEnt}) fitDV(*d);
  drawDV(dP, 400, 40);
  drawDV(dEta, 250, 40);
  drawDV(dPhi, 250, 100);
  drawDV(dX, 250, 40);
  drawDV(dVz, 250, 40);
  drawDV(dNt, 250, 40);
  drawDV(dEnt, 250, 40);
  auto dnum = [&](DV &d, const std::string &tag) {
    double mn, mx;
    for (int a = 0; a < 2; ++a)
    {
      std::string A = a ? "Z" : "Xy";
      minmax(d.gS[a], mn, mx);
      numf("eSig" + A + tag + "Min", mn, 0);
      numf("eSig" + A + tag + "Max", mx, 0);
      minmax(d.gM[a], mn, mx);
      numf("eMean" + A + tag + "Min", mn, 0);
      numf("eMean" + A + tag + "Max", mx, 0);
    }
  };
  dnum(dP, "P");
  dnum(dEta, "Eta");
  dnum(dPhi, "Phi");
  dnum(dX, "Xing");
  dnum(dVz, "Vz");
  dnum(dNt, "Ntrk");
  dnum(dEnt, "Entry");
  {
    TF1 lin("lin", "pol1", -15, 15);
    dVz.gM[0]->Fit(&lin, "QNR0");
    numf("eSlopeXyVz", lin.GetParameter(1), 2);
    numf("eESlopeXyVz", lin.GetParError(1), 2);
  }
  numf("eSigXyPfirst", dP.gS[0]->GetY()[0], 0);
  numf("eSigZPfirst", dP.gS[1]->GetY()[0], 0);
  numf("eSigXyPlast", dP.gS[0]->GetY()[dP.gS[0]->GetN() - 1], 0);
  numf("eSigZPlast", dP.gS[1]->GetY()[dP.gS[1]->GetN() - 1], 0);
  // crossing 0 vs others (bin containing 0 vs average of others), sigma xy
  {
    int b0 = dX.h[0]->GetXaxis()->FindBin(0.);
    TH1D *p0 = dX.h[0]->ProjectionY("px0", b0, b0);
    TH1D *pa = dX.h[0]->ProjectionY("pxa", 1, dX.h[0]->GetNbinsX());
    pa->Add(p0, -1);
    numf("eSigXyXbinZero", core(p0).sig * 1e4, 0);
    numf("eSigXyXbinOther", core(pa).sig * 1e4, 0);
  }

  num("eNChecks", std::to_string(nCheck));
  num("eNFailed", std::to_string(nFail));
  tex.close();
  txt.close();
  {
    std::ifstream in(note + "/tables/effdca_numbers.tex");
    std::ofstream out(note + "/tables/effdca_defaults.tex");
    out << "% generated by scripts/qa_effdca.C: fallbacks for tables/effdca_numbers.tex\n";
    std::string l;
    std::smatch m;
    while (std::getline(in, l))
      if (std::regex_search(l, m, std::regex("\\\\renewcommand\\{\\\\([A-Za-z]+)\\}"))) out << "\\newcommand{\\" << m[1] << "}{??}\n";
  }
  std::cout << "qa_effdca: failed checks " << nFail << " of " << nCheck << std::endl;
}
