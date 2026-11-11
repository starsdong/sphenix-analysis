// make_qa.C : figures, numbers and consistency checks for the TrackQA note.
// Run locally (user's ROOT style is loaded by ~/.rootlogon.C):
//   root -l -b -q 'scripts/make_qa.C("<note root>")'
// Inputs (paths relative to the trackqa project root = <note root>/..):
//   ntuples/TrackQA_00079509_00000.root            ntuple from the production DST
//   ntuples/fromtrim/TrackQA_00079509_00000.root   ntuple re-made from the trimmed DST
//   ../ntuples/KFP_EventInfo_run79509_seg000-030.root   EventInfo tree of the KFParticle test (independent code)
//   <note>/logs/dataset_software.txt, <note>/logs/run_trackqa_79509_00000.log
// Outputs: <note>/figures/*.pdf|png, <note>/tables/qa_numbers.tex, <note>/tables/tab_*.tex, <note>/logs/qa_numbers.txt
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
#include <TStyle.h>
#include <TSystem.h>
#include <TROOT.h>
#include <TBranch.h>
#include <array>
#include <cmath>
#include <TMath.h>
#include <fstream>
#include <iostream>
#include <map>
#include <regex>
#include <sstream>
#include <string>
#include <vector>

namespace qa
{
  std::string note, projDir;
  std::ofstream tex, txt;
  int nFail = 0, nCheck = 0;
  const char *sampleTag = "run3pp 79509-00000, 989 events";

  // letters-only macro names
  void num(const std::string &name, const std::string &val)
  {
    tex << "\\renewcommand{\\" << name << "}{" << val << "}\n";
    txt << name << " = " << val << "\n";
  }
  void numf(const std::string &name, double v, int prec)
  {
    std::ostringstream s;
    s.setf(std::ios::fixed);
    s.precision(prec);
    s << v;
    num(name, s.str());
  }
  void check(const std::string &name, bool ok, const std::string &what)
  {
    num(name, ok ? "PASS" : "FAIL");
    txt << "  CHECK " << name << ": " << what << " -> " << (ok ? "PASS" : "FAIL") << "\n";
    if (!ok) ++nFail;
    ++nCheck;
  }
  // defaults so every \renewcommand has a \newcommand in the preamble (written to tables/qa_defaults.tex)
  std::vector<std::string> names;

  void save(TCanvas *c, const std::string &n)
  {
    c->SaveAs((note + "/figures/" + n + ".pdf").c_str());
    c->SaveAs((note + "/figures/" + n + ".png").c_str());
  }
  void tagline(const char *extra = "")
  {
    TLatex l;
    l.SetNDC();
    l.SetTextFont(42);
    l.SetTextSize(0.042);
    l.DrawLatex(0.18, 0.935, Form("%s%s", sampleTag, extra));
  }
  TCanvas *canvas(const char *n, int w = 800, int h = 600, bool right = false)
  {
    auto *c = new TCanvas(n, n, w, h);
    c->SetTopMargin(0.09);
    c->SetRightMargin(right ? 0.15 : 0.04);
    return c;
  }
  void styleH(TH1 *h, int col, int ms = 20, int ls = 1)
  {
    h->SetLineColor(col);
    h->SetMarkerColor(col);
    h->SetMarkerStyle(ms);
    h->SetMarkerSize(1.0);
    h->SetLineStyle(ls);
    h->SetLineWidth(2);
  }
}  // namespace qa
using namespace qa;

// ---- track ntuple reader
struct Trk
{
  Int_t run, event, id;
  Long64_t ievt;
  Float_t pt, eta, phi, pcax, pcay, pcaz, chi2, quality, tvx, tvy, tvz, dcaxy, dcaxyErr, dcaz, dcazErr;
  Short_t charge, ndf, crossing, vtxMatch, tvntrk;
  UChar_t nmvtx, nintt, ntpc, ntpot;
  Bool_t hasSiSeed, hasTpcSeed;
  Int_t vtxId;
  void attach(TTree *t)
  {
    t->SetBranchAddress("run", &run);
    t->SetBranchAddress("event", &event);
    t->SetBranchAddress("ievt", &ievt);
    t->SetBranchAddress("id", &id);
    t->SetBranchAddress("pt", &pt);
    t->SetBranchAddress("eta", &eta);
    t->SetBranchAddress("phi", &phi);
    t->SetBranchAddress("charge", &charge);
    t->SetBranchAddress("pcax", &pcax);
    t->SetBranchAddress("pcay", &pcay);
    t->SetBranchAddress("pcaz", &pcaz);
    t->SetBranchAddress("chi2", &chi2);
    t->SetBranchAddress("ndf", &ndf);
    t->SetBranchAddress("quality", &quality);
    t->SetBranchAddress("crossing", &crossing);
    t->SetBranchAddress("nmvtx", &nmvtx);
    t->SetBranchAddress("nintt", &nintt);
    t->SetBranchAddress("ntpc", &ntpc);
    t->SetBranchAddress("ntpot", &ntpot);
    t->SetBranchAddress("hasSiSeed", &hasSiSeed);
    t->SetBranchAddress("hasTpcSeed", &hasTpcSeed);
    t->SetBranchAddress("vtxMatch", &vtxMatch);
    t->SetBranchAddress("vtxId", &vtxId);
    t->SetBranchAddress("tvx", &tvx);
    t->SetBranchAddress("tvy", &tvy);
    t->SetBranchAddress("tvz", &tvz);
    t->SetBranchAddress("tvntrk", &tvntrk);
    t->SetBranchAddress("dcaxy", &dcaxy);
    t->SetBranchAddress("dcaxyErr", &dcaxyErr);
    t->SetBranchAddress("dcaz", &dcaz);
    t->SetBranchAddress("dcazErr", &dcazErr);
  }
};
struct Evt
{
  Int_t run, event, nTrk, nTrkSil, nTrkINTT, nTrkMVTXonly, nTrkTPConly, nTrkX0, nTrkSilX0, nXingTrk, nVtx, nVtxX0, vntrk, GL1_bunch;
  Long64_t GL1_BCO;
  ULong64_t GL1_scaled, GL1_live;
  Float_t vx, vy, vz, vchi2ndf;
  void attach(TTree *t)
  {
    t->SetBranchAddress("run", &run);
    t->SetBranchAddress("event", &event);
    t->SetBranchAddress("GL1_BCO", &GL1_BCO);
    t->SetBranchAddress("GL1_scaled", &GL1_scaled);
    t->SetBranchAddress("GL1_live", &GL1_live);
    t->SetBranchAddress("GL1_bunch", &GL1_bunch);
    t->SetBranchAddress("nTrk", &nTrk);
    t->SetBranchAddress("nTrkSil", &nTrkSil);
    t->SetBranchAddress("nTrkINTT", &nTrkINTT);
    t->SetBranchAddress("nTrkMVTXonly", &nTrkMVTXonly);
    t->SetBranchAddress("nTrkTPConly", &nTrkTPConly);
    t->SetBranchAddress("nTrkX0", &nTrkX0);
    t->SetBranchAddress("nTrkSilX0", &nTrkSilX0);
    t->SetBranchAddress("nXingTrk", &nXingTrk);
    t->SetBranchAddress("nVtx", &nVtx);
    t->SetBranchAddress("nVtxX0", &nVtxX0);
    t->SetBranchAddress("vx", &vx);
    t->SetBranchAddress("vy", &vy);
    t->SetBranchAddress("vz", &vz);
    t->SetBranchAddress("vntrk", &vntrk);
    t->SetBranchAddress("vchi2ndf", &vchi2ndf);
  }
};

std::string fmtInt(long v)
{
  std::string s = std::to_string(v), o;
  int n = 0;
  for (int i = s.size() - 1; i >= 0; --i)
  {
    o.insert(o.begin(), s[i]);
    if (++n % 3 == 0 && i > 0 && s[i - 1] != '-') o.insert(o.begin(), ',');
  }
  return o;
}

// Gaussian core width: iterative fit in +-2 sigma
std::pair<double, double> coreSigma(TH1 *h)
{
  if (h->GetEntries() < 50) return {-1, -1};
  double m = h->GetBinCenter(h->GetMaximumBin()), s = h->GetRMS();
  TF1 f("g", "gaus", m - 2 * s, m + 2 * s);
  for (int it = 0; it < 4; ++it)
  {
    f.SetRange(m - 2 * s, m + 2 * s);
    f.SetParameters(h->GetMaximum(), m, s);
    h->Fit(&f, "QNR0");
    m = f.GetParameter(1);
    s = std::fabs(f.GetParameter(2));
  }
  return {s, f.GetParError(2)};
}

void make_qa(const char *noteRoot = ".")
{
  note = noteRoot;
  projDir = note + "/..";
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  TH1::SetDefaultSumw2();
  gSystem->mkdir((note + "/figures").c_str(), true);
  gSystem->mkdir((note + "/tables").c_str(), true);
  tex.open(note + "/tables/qa_numbers.tex");
  txt.open(note + "/logs/qa_numbers.txt");
  tex << "% generated by scripts/make_qa.C -- do not edit\n";

  // ------------------------------------------------------------------ inputs
  TFile *f = TFile::Open((projDir + "/ntuples/TrackQA_00079509_00000.root").c_str());
  TFile *ft = TFile::Open((projDir + "/ntuples/fromtrim/TrackQA_00079509_00000.root").c_str());
  TFile *fk = TFile::Open((projDir + "/../ntuples/KFP_EventInfo_run79509_seg000-030.root").c_str());
  TTree *te = (TTree *) f->Get("event"), *tt = (TTree *) f->Get("track");
  TTree *te2 = (TTree *) ft->Get("event"), *tt2 = (TTree *) ft->Get("track");
  TTree *tk = (TTree *) fk->Get("EventInfo");
  const long nEv = te->GetEntries(), nTr = tt->GetEntries();

  // dataset/software log
  long dstEntries = -1, dstBytes = -1;
  std::map<std::string, long> dstBr, trimBr;
  {
    std::ifstream in(note + "/logs/dataset_software.txt");
    std::string l;
    while (std::getline(in, l))
    {
      std::istringstream is(l);
      std::string a, b;
      long v;
      is >> a;
      if (a == "DST_ENTRIES") is >> dstEntries;
      if (a == "DSTBR" && (is >> b >> v)) dstBr[b] = v;
      if (a == "TRIMBR" && (is >> b >> v)) trimBr[b] = v;
      if (a.size() > 3 && a[0] == '-' && l.find("DST_TRKR_TRACKS_run3pp") != std::string::npos)
      {
        std::istringstream js(l);
        std::string x;
        for (int k = 0; k < 5; ++k) js >> x;
        dstBytes = std::stol(x);
      }
    }
  }
  // job log: cpu / wall / rss / module time / vtxMatch counts
  double userT = -1, sysT = -1, modMs = -1;
  std::string wall;
  long rss = -1, logMatch[4] = {-1, -1, -1, -1}, logEv = -1, logTr = -1;
  {
    std::ifstream in(note + "/logs/run_trackqa_79509_00000.log");
    std::string l;
    std::smatch m;
    while (std::getline(in, l))
    {
      if (std::regex_search(l, m, std::regex("User time \\(seconds\\): ([0-9.]+)"))) userT = std::stod(m[1]);
      if (std::regex_search(l, m, std::regex("System time \\(seconds\\): ([0-9.]+)"))) sysT = std::stod(m[1]);
      if (std::regex_search(l, m, std::regex("Elapsed .*: ([0-9:.]+)"))) wall = m[1];
      if (std::regex_search(l, m, std::regex("Maximum resident set size \\(kbytes\\): ([0-9]+)"))) rss = std::stol(m[1]);
      if (std::regex_search(l, m, std::regex("TrackQATree_TOP: per event time \\(ms\\): +([0-9.]+)"))) modMs = std::stod(m[1]);
      if (std::regex_search(l, m, std::regex("TrackQATree: events ([0-9]+) tracks written ([0-9]+) vtxMatch\\(0,1,2,3\\) ([0-9]+) ([0-9]+) ([0-9]+) ([0-9]+)")))
      {
        logEv = std::stol(m[1]);
        logTr = std::stol(m[2]);
        for (int k = 0; k < 4; ++k) logMatch[k] = std::stol(m[3 + k]);
      }
    }
  }
  long fsNt = f->GetSize(), fsTrim = 0;
  {
    TFile *fd = TFile::Open((projDir + "/ntuples/DST_TRKR_TRACKS_TRIM-00079509-00000.root").c_str());
    fsTrim = fd->GetSize();
    fd->Close();
  }

  // ------------------------------------------------------------------ event loop
  Evt e;
  e.attach(te);
  TH1D *hNtrk = new TH1D("hNtrk", ";tracks per event;events", 60, 0, 300);
  TH1D *hNsil = new TH1D("hNsil", "", 60, 0, 300);
  TH1D *hNtpc = new TH1D("hNtpc", "", 60, 0, 300);
  TH1D *hNx0 = new TH1D("hNx0", "", 60, 0, 300);
  TH1D *hNvtx = new TH1D("hNvtx", ";vertices per event;events", 40, -0.5, 39.5);
  TH1D *hNvtx0 = new TH1D("hNvtx0", "", 40, -0.5, 39.5);
  TH1D *hNxing = new TH1D("hNxing", ";crossings with silicon tracks;events", 50, -0.5, 49.5);
  TH1D *hVz = new TH1D("hVz", ";#it{v}_{z} of crossing-0 vertex (cm);events", 50, -50, 50);
  TH2D *hVxy = new TH2D("hVxy", ";#it{v}_{x} (cm);#it{v}_{y} (cm)", 60, -0.3, 0.3, 60, -0.3, 0.3);
  TH1D *hVntrk = new TH1D("hVntrk", ";tracks in crossing-0 vertex;events", 40, 0.5, 40.5);
  long sumTrk = 0, sumSil = 0, sumINTT = 0, sumMonly = 0, sumTpc = 0, sumX0 = 0, sumSilX0 = 0, sumVtx = 0, sumVtx0 = 0, nEvVtx0 = 0, sumXing = 0;
  bool typeSumOK = true, silSumOK = true;
  std::vector<int> evNumbers(nEv);
  std::vector<std::array<int, 6>> evCounts(nEv);
  int run0 = -1;
  for (long i = 0; i < nEv; ++i)
  {
    te->GetEntry(i);
    run0 = e.run;
    evNumbers[i] = e.event;
    evCounts[i] = {e.nTrk, e.nTrkINTT, e.nTrkMVTXonly, e.nTrkTPConly, e.nTrkX0, e.nVtx};
    sumTrk += e.nTrk;
    sumSil += e.nTrkSil;
    sumINTT += e.nTrkINTT;
    sumMonly += e.nTrkMVTXonly;
    sumTpc += e.nTrkTPConly;
    sumX0 += e.nTrkX0;
    sumSilX0 += e.nTrkSilX0;
    sumVtx += e.nVtx;
    sumVtx0 += e.nVtxX0;
    sumXing += e.nXingTrk;
    if (e.nTrkINTT + e.nTrkMVTXonly + e.nTrkTPConly != e.nTrk) typeSumOK = false;
    if (e.nTrkINTT + e.nTrkMVTXonly != e.nTrkSil) silSumOK = false;
    hNtrk->Fill(e.nTrk);
    hNsil->Fill(e.nTrkSil);
    hNtpc->Fill(e.nTrkTPConly);
    hNx0->Fill(e.nTrkSilX0);
    hNvtx->Fill(e.nVtx);
    hNvtx0->Fill(e.nVtxX0);
    hNxing->Fill(e.nXingTrk);
    if (e.nVtxX0 > 0)
    {
      ++nEvVtx0;
      hVz->Fill(e.vz);
      hVxy->Fill(e.vx, e.vy);
      hVntrk->Fill(e.vntrk);
    }
  }

  // ------------------------------------------------------------------ track loop
  Trk t;
  t.attach(tt);
  TH1D *hPt[3], *hEta[3], *hPhi[3];
  const char *cls[3] = {"INTT", "MVTXonly", "TPConly"};
  for (int k = 0; k < 3; ++k)
  {
    hPt[k] = new TH1D(Form("hPt%d", k), ";#it{p}_{T} (GeV/#it{c});tracks / event / (GeV/#it{c})", 100, 0, 5);
    hEta[k] = new TH1D(Form("hEta%d", k), ";#eta;tracks / event / 0.05", 80, -2, 2);
    hPhi[k] = new TH1D(Form("hPhi%d", k), ";#phi (rad);tracks / event / rad", 72, -TMath::Pi(), TMath::Pi());
  }
  TH2D *hEtaPhi = new TH2D("hEtaPhi", ";#eta;#phi (rad)", 60, -1.5, 1.5, 72, -TMath::Pi(), TMath::Pi());
  TH2D *hEtaPhiT = new TH2D("hEtaPhiT", ";#eta;#phi (rad)", 60, -1.5, 1.5, 72, -TMath::Pi(), TMath::Pi());
  TH1D *hNm = new TH1D("hNm", ";clusters on track;tracks", 61, -0.5, 60.5);
  TH1D *hNi = new TH1D("hNi", "", 61, -0.5, 60.5);
  TH1D *hNt = new TH1D("hNt", "", 61, -0.5, 60.5);
  TH1D *hNp = new TH1D("hNp", "", 61, -0.5, 60.5);
  TH1D *hQ[2];
  hQ[0] = new TH1D("hQ0", ";#chi^{2}/ndf;tracks (normalised)", 100, 0, 20);
  hQ[1] = new TH1D("hQ1", "", 100, 0, 20);
  TH1D *hX = new TH1D("hX", ";beam crossing of silicon track;tracks / event", 130, -150, 500);
  TH1D *hDxy[4], *hDz[4];
  for (int k = 0; k < 4; ++k)
  {
    hDxy[k] = new TH1D(Form("hDxy%d", k), ";DCA_{xy} (cm);tracks (normalised)", 100, -1, 1);
    hDz[k] = new TH1D(Form("hDz%d", k), ";DCA_{z} (cm);tracks (normalised)", 100, -2, 2);
  }
  const int nPtB = 8;
  double ptB[nPtB + 1] = {0.2, 0.3, 0.4, 0.5, 0.7, 1.0, 1.5, 2.0, 3.0};
  TH1D *hRes[2][2][nPtB];  // [xy,z][vtxMatch1, vtxMatch1 & vertex ntrk>=10][pt]
  for (int a = 0; a < 2; ++a)
    for (int b = 0; b < 2; ++b)
      for (int k = 0; k < nPtB; ++k) hRes[a][b][k] = new TH1D(Form("hRes%d%d%d", a, b, k), "", 400, -0.2, 0.2);
  long nTpcOnlyNonX0 = 0, nTpcOnly = 0, nNoSeedSi = 0, nSil = 0, nGood = 0, nMatch[4] = {0, 0, 0, 0};
  long nDcaFinite = 0, nWithVtx = 0, nMatchedTpcOnly = 0, nPosQ = 0, nNegQ = 0;
  long nNdfNeg = 0, nTpcMatch12 = 0, nSilMatch3 = 0;
  int ntpcMin = 999, ntpcMax = -1, nmvtxMax = -1, ninttMax = -1, ntpotMax = -1;
  long nTpot = 0, nPt02 = 0;
  for (long i = 0; i < nTr; ++i)
  {
    tt->GetEntry(i);
    const bool sil = (t.nmvtx + t.nintt) > 0;
    const int c = t.nintt > 0 ? 0 : (t.nmvtx > 0 ? 1 : 2);
    if (t.vtxMatch >= 0 && t.vtxMatch < 4) ++nMatch[t.vtxMatch];
    if (t.charge > 0) ++nPosQ; else ++nNegQ;
    hPt[c]->Fill(t.pt);
    if (t.pt > 0.2)
    {
      hEta[c]->Fill(t.eta);
      hPhi[c]->Fill(t.phi);
    }
    if (sil)
    {
      ++nSil;
      if (!t.hasSiSeed) ++nNoSeedSi;
      if (t.vtxMatch == 3) ++nSilMatch3;
      hX->Fill(t.crossing);
      if (t.pt > 0.2) hEtaPhi->Fill(t.eta, t.phi);
    }
    else
    {
      ++nTpcOnly;
      if (t.crossing != 0) ++nTpcOnlyNonX0;
      if (t.vtxMatch > 0) ++nMatchedTpcOnly;
      if (t.vtxMatch == 1 || t.vtxMatch == 2) ++nTpcMatch12;
      if (t.pt > 0.2) hEtaPhiT->Fill(t.eta, t.phi);
    }
    ntpcMin = std::min<int>(ntpcMin, t.ntpc);
    ntpcMax = std::max<int>(ntpcMax, t.ntpc);
    nmvtxMax = std::max<int>(nmvtxMax, t.nmvtx);
    ninttMax = std::max<int>(ninttMax, t.nintt);
    ntpotMax = std::max<int>(ntpotMax, t.ntpot);
    if (t.ntpot > 0) ++nTpot;
    if (t.pt < 0.2) ++nPt02;
    hNm->Fill(t.nmvtx);
    hNi->Fill(t.nintt);
    hNt->Fill(t.ntpc);
    hNp->Fill(t.ntpot);
    if (t.ndf > 0) hQ[sil ? 0 : 1]->Fill(t.chi2 / t.ndf);
    else ++nNdfNeg;
    if (t.vtxMatch > 0)
    {
      ++nWithVtx;
      if (std::isfinite(t.dcaxy) && std::isfinite(t.dcaz)) ++nDcaFinite;
      if (t.pt > 0.5)  // vtxMatch 1,2 are silicon tracks, 3 TPC-only (check chkTpcVtxId)
      {
        hDxy[t.vtxMatch]->Fill(t.dcaxy);
        hDz[t.vtxMatch]->Fill(t.dcaz);
      }
    }
    // "good" track: MVTX>=3, INTT>=1, TPC>=20 (DCA resolution sample)
    const bool good = t.nmvtx >= 3 && t.nintt >= 1 && t.ntpc >= 20;
    if (good) ++nGood;
    if (good && t.vtxMatch == 1)
    {
      for (int k = 0; k < nPtB; ++k)
        if (t.pt >= ptB[k] && t.pt < ptB[k + 1])
        {
          hRes[0][0][k]->Fill(t.dcaxy);
          hRes[1][0][k]->Fill(t.dcaz);
          if (t.tvntrk >= 10)
          {
            hRes[0][1][k]->Fill(t.dcaxy);
            hRes[1][1][k]->Fill(t.dcaz);
          }
        }
    }
  }

  // ------------------------------------------------------------------ consistency checks
  // (1) trimmed-DST ntuple identical to the production-DST ntuple
  bool sameTrim = (te2->GetEntries() == nEv) && (tt2->GetEntries() == nTr);
  long nDiff = 0;
  if (sameTrim)
  {
    Trk t2;
    t2.attach(tt2);
    for (long i = 0; i < nTr; ++i)
    {
      tt->GetEntry(i);
      tt2->GetEntry(i);
      auto eqf = [](float a, float b) { return a == b || (std::isnan(a) && std::isnan(b)); };
      if (t.event != t2.event || t.id != t2.id || !eqf(t.pt, t2.pt) || !eqf(t.eta, t2.eta) || !eqf(t.phi, t2.phi) ||
          !eqf(t.dcaxy, t2.dcaxy) || !eqf(t.dcaz, t2.dcaz) || t.crossing != t2.crossing || t.nmvtx != t2.nmvtx ||
          t.nintt != t2.nintt || t.ntpc != t2.ntpc || t.vtxMatch != t2.vtxMatch || !eqf(t.chi2, t2.chi2))
        ++nDiff;
    }
    Evt e2;
    e2.attach(te2);
    for (long i = 0; i < nEv; ++i)
    {
      te->GetEntry(i);
      te2->GetEntry(i);
      if (e.event != e2.event || e.nTrk != e2.nTrk || e.nVtx != e2.nVtx || e.GL1_BCO != e2.GL1_BCO || e.vz != e2.vz) ++nDiff;
    }
  }
  // (2) per-event counts vs KFParticle-test EventInfo (independent module, same DST)
  std::map<int, long> kIdx;
  Int_t kRun, kEvt, kN[6];
  tk->SetBranchAddress("runNumber", &kRun);
  tk->SetBranchAddress("eventNumber", &kEvt);
  tk->SetBranchAddress("nTrk", &kN[0]);
  tk->SetBranchAddress("nTrkINTT", &kN[1]);
  tk->SetBranchAddress("nTrkMVTXonly", &kN[2]);
  tk->SetBranchAddress("nTrkTPConly", &kN[3]);
  tk->SetBranchAddress("nTrkX0", &kN[4]);
  tk->SetBranchAddress("nVtx", &kN[5]);
  for (long i = 0; i < tk->GetEntries(); ++i)
  {
    tk->GetEntry(i);
    if (kRun == run0 && !kIdx.count(kEvt)) kIdx[kEvt] = i;  // first occurrence
  }
  long kFound = 0, kAgree = 0;
  for (long i = 0; i < nEv; ++i)
  {
    auto it = kIdx.find(evNumbers[i]);
    if (it == kIdx.end()) continue;
    ++kFound;
    tk->GetEntry(it->second);
    bool ok = true;
    for (int k = 0; k < 6; ++k)
      if (kN[k] != evCounts[i][k]) ok = false;
    if (ok) ++kAgree;
  }

  // ------------------------------------------------------------------ numbers
  for (const char *n : {"nEvents", "nTracks"}) (void) n;
  num("runNo", std::to_string(run0));
  num("nEvents", fmtInt(nEv));
  num("nTracks", fmtInt(nTr));
  num("nDSTEntries", fmtInt(dstEntries));
  numf("trkPerEv", double(sumTrk) / nEv, 1);
  numf("silPerEv", double(sumSil) / nEv, 1);
  numf("inttPerEv", double(sumINTT) / nEv, 1);
  numf("mvtxOnlyPerEv", double(sumMonly) / nEv, 1);
  numf("tpcOnlyPerEv", double(sumTpc) / nEv, 1);
  numf("silXzeroPerEv", double(sumSilX0) / nEv, 2);
  numf("trkXzeroPerEv", double(sumX0) / nEv, 1);
  numf("fracTpcOnly", 100. * sumTpc / sumTrk, 1);
  numf("fracINTT", 100. * sumINTT / sumTrk, 1);
  numf("fracMVTXonly", 100. * sumMonly / sumTrk, 1);
  numf("vtxPerEv", double(sumVtx) / nEv, 2);
  numf("vtxXzeroPerEv", double(sumVtx0) / nEv, 3);
  numf("fracEvVtxXzero", 100. * nEvVtx0 / nEv, 1);
  numf("xingPerEv", double(sumXing) / nEv, 1);
  numf("vzMean", hVz->GetMean(), 2);
  numf("vzRMS", hVz->GetRMS(), 1);
  numf("vxMean", hVxy->GetMean(1) * 1e4, 0);
  numf("vyMean", hVxy->GetMean(2) * 1e4, 0);
  numf("vntrkMean", hVntrk->GetMean(), 1);
  numf("fracPos", 100. * nPosQ / nTr, 1);
  for (int k = 0; k < 4; ++k)
  {
    const char *w[4] = {"Zero", "One", "Two", "Three"};
    num(std::string("nMatch") + w[k], fmtInt(nMatch[k]));
    numf(std::string("fracMatch") + w[k], 100. * nMatch[k] / nTr, 1);
  }
  numf("fracGood", 100. * nGood / nTr, 1);
  numf("sizeDSTkB", dstBytes / 1024. / nEv, 1);
  numf("sizeTrimkB", fsTrim / 1024. / nEv, 1);
  numf("sizeNtkB", fsNt / 1024. / nEv, 2);
  numf("sizeNtMB", fsNt / 1024. / 1024., 2);
  numf("sizeTrimMB", fsTrim / 1024. / 1024., 1);
  numf("sizeDSTMB", dstBytes / 1024. / 1024., 1);
  numf("ratioNt", double(dstBytes) / fsNt, 0);
  numf("ratioTrim", double(dstBytes) / fsTrim, 1);
  // projection to the full run3pp sample (~5e8 events, from the KFParticle note / README)
  const double nFull = 5e8, TB = 1024. * 1024. * 1024. * 1024.;
  numf("projNtFullTB", fsNt / double(nEv) * nFull / TB, 1);
  numf("projNtTenPctTB", fsNt / double(nEv) * nFull / TB / 10., 2);
  numf("projTrimFullTB", fsTrim / double(nEv) * nFull / TB, 1);
  numf("ratioTrimNt", double(fsTrim) / fsNt, 1);
  numf("bytesPerTrk", double(tt->GetZipBytes()) / nTr, 1);
  numf("cpuUser", userT, 1);
  numf("cpuSys", sysT, 1);
  num("wallTime", wall);
  numf("rssMB", rss / 1024., 0);
  numf("modMsPerEv", modMs, 2);
  numf("cpuMsPerEv", 1000. * (userT + sysT) / nEv, 1);
  // cluster/seed related DST fractions
  long dstTot = 0;
  for (auto &[k, v] : dstBr) dstTot += v;
  auto pct = [&](const std::string &b) { return 100. * dstBr[b] / dstTot; };
  numf("pctClus", pct("DST#TRKR#TRKR_CLUSTER_SEED"), 1);
  numf("pctSeeds", pct("DST#SVTX#SiliconTrackSeedContainer") + pct("DST#SVTX#TpcTrackSeedContainer") + pct("DST#SVTX#SvtxTrackSeedContainer"), 1);
  numf("pctTrkMap", pct("DST#SVTX#SvtxTrackMap"), 1);
  numf("pctVtxMap", pct("DST#SVTX#SvtxVertexMap"), 2);
  numf("trkMapkBPerEv", dstBr["DST#SVTX#SvtxTrackMap"] / 1024. / nEv, 1);

  // checks
  check("chkEvEntries", nEv == dstEntries && logEv == nEv, "event-tree entries == DST entries == events in job log");
  check("chkTrkEntries", sumTrk == nTr && logTr == nTr, "sum of nTrk over events == track-tree entries == job log (minPt=0)");
  check("chkTypeSum", typeSumOK && silSumOK, "nTrkINTT+nTrkMVTXonly+nTrkTPConly == nTrk and nTrkINTT+nTrkMVTXonly == nTrkSil in every event");
  check("chkTpcOnlyXzero", nTpcOnlyNonX0 == 0, "all TPC-only tracks have crossing 0");
  check("chkMatchLog", logMatch[0] == nMatch[0] && logMatch[1] == nMatch[1] && logMatch[2] == nMatch[2] && logMatch[3] == nMatch[3], "vtxMatch counts in ntuple == job log");
  check("chkDcaFinite", nDcaFinite == nWithVtx, "DCA finite for every track with an associated vertex");
  check("chkSiSeed", nNoSeedSi == 0, "every silicon track has a silicon seed (id)");
  check("chkTrim", sameTrim && nDiff == 0, "ntuple from trimmed DST identical to ntuple from production DST (entry by entry)");
  check("chkKFP", kFound == nEv && kAgree == nEv, "per-event nTrk, nTrkINTT, nTrkMVTXonly, nTrkTPConly, nTrkX0, nVtx == KFP-test EventInfo");
  num("nKFPFound", fmtInt(kFound));
  num("ntpcMin", std::to_string(ntpcMin));
  num("ntpcMax", std::to_string(ntpcMax));
  num("nmvtxMax", std::to_string(nmvtxMax));
  num("ninttMax", std::to_string(ninttMax));
  num("ntpotMax", std::to_string(ntpotMax));
  numf("fracTpot", 100. * nTpot / nTr, 1);
  numf("fracPtLow", 100. * nPt02 / nTr, 1);
  numf("ntpcMean", hNt->GetMean(), 1);
  numf("qualSilMean", hQ[0]->GetMean(), 2);
  numf("qualTpcMean", hQ[1]->GetMean(), 2);
  num("nKFPAgree", fmtInt(kAgree));
  num("nTrimDiff", std::to_string(nDiff));
  num("nNdfZero", std::to_string(nNdfNeg));
  num("nMatchedTpcOnly", fmtInt(nMatchedTpcOnly));
  check("chkTpcVtxId", nTpcMatch12 == 0 && nSilMatch3 == 0, "TPC-only tracks never carry a vertex id (no vtxMatch 1/2) and silicon tracks never need the nearest-vertex fallback (no vtxMatch 3)");
  numf("fracMatchedTpcOnly", 100. * nMatchedTpcOnly / std::max(1L, nTpcOnly), 1);

  // ------------------------------------------------------------------ figures
  // event multiplicities
  {
    auto *c = canvas("c1");
    c->SetLogy();
    styleH(hNtrk, kBlack);
    styleH(hNsil, kBlue + 1);
    styleH(hNtpc, kRed + 1);
    styleH(hNx0, kGreen + 2);
    hNtrk->SetMinimum(0.5);
    hNtrk->SetMaximum(hNtrk->GetMaximum() * 20);
    hNtrk->Draw("hist");
    hNsil->Draw("hist same");
    hNtpc->Draw("hist same");
    hNx0->Draw("hist same");
    auto *lg = new TLegend(0.55, 0.66, 0.95, 0.89);
    lg->SetBorderSize(0);
    lg->SetTextSize(0.04);
    lg->AddEntry(hNtrk, "all tracks", "l");
    lg->AddEntry(hNsil, "with silicon (MVTX/INTT)", "l");
    lg->AddEntry(hNtpc, "TPC-only", "l");
    lg->AddEntry(hNx0, "silicon, crossing 0", "l");
    lg->Draw();
    tagline();
    save(c, "fig_evt_ntrk");
  }
  {
    auto *c = canvas("c2", 1200, 500);
    c->Divide(2, 1);
    c->cd(1);
    gPad->SetTopMargin(0.09);
    gPad->SetRightMargin(0.04);
    gPad->SetLogy();
    styleH(hNvtx, kBlack);
    styleH(hNvtx0, kRed + 1);
    hNvtx->SetMinimum(0.5);
    hNvtx->SetMaximum(hNvtx0->GetMaximum() * 5);
    hNvtx->Draw("hist");
    hNvtx0->Draw("hist same");
    auto *lg = new TLegend(0.45, 0.72, 0.95, 0.89);
    lg->SetBorderSize(0);
    lg->SetTextSize(0.045);
    lg->AddEntry(hNvtx, "all vertices", "l");
    lg->AddEntry(hNvtx0, "crossing-0 vertices", "l");
    lg->Draw();
    tagline();
    c->cd(2);
    gPad->SetTopMargin(0.09);
    gPad->SetRightMargin(0.04);
    styleH(hVz, kBlack);
    hVz->Draw("e");
    save(c, "fig_evt_vtx");
  }
  {
    auto *c = canvas("c3", 1200, 500);
    c->Divide(2, 1);
    c->cd(1);
    gPad->SetTopMargin(0.09);
    gPad->SetRightMargin(0.04);
    hX->Scale(1. / nEv);
    styleH(hX, kBlue + 1);
    hX->Draw("hist");
    tagline();
    c->cd(2);
    gPad->SetTopMargin(0.09);
    gPad->SetRightMargin(0.04);
    styleH(hNxing, kBlack);
    hNxing->Draw("hist");
    save(c, "fig_evt_crossing");
  }
  // track kinematics
  {
    auto *c = canvas("c4", 1500, 540);
    c->Divide(3, 1);
    int col[3] = {kBlue + 1, kGreen + 2, kRed + 1};
    const char *lab[3] = {"INTT (#geq1 INTT cluster)", "MVTX, no INTT", "TPC-only"};
    TH1D **hs[3] = {hPt, hEta, hPhi};
    for (int p = 0; p < 3; ++p)
    {
      c->cd(p + 1);
      gPad->SetTopMargin(0.12);
      gPad->SetRightMargin(0.04);
      gPad->SetLeftMargin(0.2);
      if (p == 0) gPad->SetLogy();
      double mx = 0;
      for (int k = 0; k < 3; ++k)
      {
        hs[p][k]->Scale(1. / nEv, "width");
        styleH(hs[p][k], col[k]);
        mx = std::max(mx, hs[p][k]->GetMaximum());
      }
      hs[p][0]->SetMaximum(p == 0 ? mx * 5 : mx * 1.45);
      if (p == 0) hs[p][0]->SetMinimum(1e-3);
      else hs[p][0]->SetMinimum(0);
      for (int k = 0; k < 3; ++k) hs[p][k]->Draw(k ? "hist same" : "hist");
      if (p == 0)
      {
        auto *lg = new TLegend(0.35, 0.7, 0.95, 0.89);
        lg->SetBorderSize(0);
        lg->SetTextSize(0.045);
        for (int k = 0; k < 3; ++k) lg->AddEntry(hs[p][k], lab[k], "l");
        lg->Draw();
        tagline();
      }
      else
      {
        TLatex l;
        l.SetNDC();
        l.SetTextSize(0.05);
        l.DrawLatex(0.25, 0.83, "#it{p}_{T} > 0.2 GeV/#it{c}");
      }
    }
    save(c, "fig_trk_kin");
  }
  {
    auto *c = canvas("c5", 1200, 500, true);
    c->Divide(2, 1);
    TH2D *hh[2] = {hEtaPhi, hEtaPhiT};
    const char *lab[2] = {"silicon tracks, #it{p}_{T} > 0.2 GeV/#it{c}", "TPC-only tracks, #it{p}_{T} > 0.2 GeV/#it{c}"};
    for (int p = 0; p < 2; ++p)
    {
      c->cd(p + 1);
      gPad->SetTopMargin(0.09);
      gPad->SetRightMargin(0.15);
      hh[p]->Draw("colz");
      TLatex l;
      l.SetNDC();
      l.SetTextSize(0.045);
      l.DrawLatex(0.18, 0.935, lab[p]);
    }
    save(c, "fig_trk_etaphi");
  }
  {
    auto *c = canvas("c6", 1200, 500);
    c->Divide(2, 1);
    c->cd(1);
    gPad->SetTopMargin(0.09);
    gPad->SetRightMargin(0.04);
    gPad->SetLogy();
    styleH(hNm, kGreen + 2);
    styleH(hNi, kBlue + 1);
    styleH(hNt, kRed + 1);
    styleH(hNp, kMagenta + 1);
    hNm->SetMinimum(10);
    hNm->SetMaximum(std::max({hNm->GetMaximum(), hNi->GetMaximum(), hNt->GetMaximum(), hNp->GetMaximum()}) * 5);
    hNm->Draw("hist");
    hNi->Draw("hist same");
    hNt->Draw("hist same");
    hNp->Draw("hist same");
    auto *lg = new TLegend(0.6, 0.66, 0.95, 0.89);
    lg->SetBorderSize(0);
    lg->SetTextSize(0.045);
    lg->AddEntry(hNm, "MVTX", "l");
    lg->AddEntry(hNi, "INTT", "l");
    lg->AddEntry(hNt, "TPC", "l");
    lg->AddEntry(hNp, "TPOT", "l");
    lg->Draw();
    tagline();
    c->cd(2);
    gPad->SetTopMargin(0.09);
    gPad->SetRightMargin(0.04);
    gPad->SetLogy();
    for (int k = 0; k < 2; ++k) hQ[k]->Scale(1. / hQ[k]->Integral());
    styleH(hQ[0], kBlue + 1);
    styleH(hQ[1], kRed + 1);
    hQ[0]->SetMaximum(std::max(hQ[0]->GetMaximum(), hQ[1]->GetMaximum()) * 5);
    hQ[0]->Draw("hist");
    hQ[1]->Draw("hist same");
    auto *lg2 = new TLegend(0.55, 0.75, 0.95, 0.89);
    lg2->SetBorderSize(0);
    lg2->SetTextSize(0.045);
    lg2->AddEntry(hQ[0], "silicon tracks", "l");
    lg2->AddEntry(hQ[1], "TPC-only tracks", "l");
    lg2->Draw();
    save(c, "fig_trk_hits");
  }
  // DCA distributions by vertex association (pT>0.5)
  {
    auto *c = canvas("c7", 1200, 500);
    c->Divide(2, 1);
    int col[4] = {kGray, kBlack, kBlue + 1, kRed + 1};
    const char *lab[4] = {"", "Si track, in vertex fit (1)", "Si track, vertex id, not in fit (2)", "TPC-only, nearest X=0 vertex (3)"};
    TH1D **hs[2] = {hDxy, hDz};
    for (int p = 0; p < 2; ++p)
    {
      c->cd(p + 1);
      gPad->SetTopMargin(0.09);
      gPad->SetRightMargin(0.04);
      gPad->SetLogy();
      double mx = 0;
      for (int k = 1; k < 4; ++k)
      {
        if (hs[p][k]->Integral() > 0) hs[p][k]->Scale(1. / hs[p][k]->Integral());
        styleH(hs[p][k], col[k]);
        mx = std::max(mx, hs[p][k]->GetMaximum());
      }
      hs[p][1]->SetMaximum(mx * 30);
      hs[p][1]->SetMinimum(1e-5);
      hs[p][1]->Draw("hist");
      for (int k : {2, 3}) hs[p][k]->Draw("hist same");
      if (p == 0)
      {
        auto *lg = new TLegend(0.2, 0.68, 0.95, 0.89);
        lg->SetBorderSize(0);
        lg->SetTextSize(0.04);
        lg->SetNColumns(1);
        for (int k : {1, 2, 3}) lg->AddEntry(hs[p][k], lab[k], "l");
        lg->Draw();
        tagline();
      }
      else
      {
        TLatex l;
        l.SetNDC();
        l.SetTextSize(0.05);
        l.DrawLatex(0.22, 0.83, "#it{p}_{T} > 0.5 GeV/#it{c}");
      }
    }
    save(c, "fig_trk_dca");
  }
  // DCA core width vs pT
  std::ofstream tab(note + "/tables/tab_dcares.tex");
  tab << "% generated by scripts/make_qa.C\n\\begin{tabular}{lrrrrrr}\n\\toprule\n"
      << "& \\multicolumn{3}{c}{all vertices} & \\multicolumn{3}{c}{vertex with $\\geq10$ tracks}\\\\\n$p_T$ (GeV/$c$) & $N$ & $\\sigma_{xy}$ ($\\mu$m) & $\\sigma_{z}$ ($\\mu$m) & $N$ & $\\sigma_{xy}$ ($\\mu$m) & $\\sigma_{z}$ ($\\mu$m)\\\\\n\\midrule\n";
  TGraphErrors *gR[2][2];
  for (int a = 0; a < 2; ++a)
    for (int b = 0; b < 2; ++b) gR[a][b] = new TGraphErrors();
  for (int k = 0; k < nPtB; ++k)
  {
    tab << Form("%.1f--%.1f", ptB[k], ptB[k + 1]);
    for (int b = 0; b < 2; ++b)
    {
      tab << " & " << fmtInt((long) hRes[0][b][k]->GetEntries());
      for (int a = 0; a < 2; ++a)
      {
        auto s = coreSigma(hRes[a][b][k]);
        if (s.first > 0)
        {
          int n = gR[a][b]->GetN();
          gR[a][b]->SetPoint(n, 0.5 * (ptB[k] + ptB[k + 1]), s.first * 1e4);
          gR[a][b]->SetPointError(n, 0, s.second * 1e4);
          tab << Form(" & %.0f $\\pm$ %.0f", s.first * 1e4, s.second * 1e4);
        }
        else
          tab << " & --";
      }
    }
    tab << "\\\\\n";
  }
  tab << "\\bottomrule\n\\end{tabular}\n";
  tab.close();
  {
    auto *c = canvas("c8");
    TH1F *fr = c->DrawFrame(0, 0, 3.2, 320);
    fr->SetXTitle("#it{p}_{T} (GeV/#it{c})");
    fr->SetYTitle("DCA Gaussian core #sigma (#mum)");
    int col[2] = {kBlack, kBlue + 1};
    int ms[2][2] = {{20, 24}, {21, 25}};
    auto *lg = new TLegend(0.4, 0.68, 0.95, 0.89);
    lg->SetBorderSize(0);
    lg->SetTextSize(0.038);
    const char *al[2] = {"DCA_{xy}", "DCA_{z}"};
    const char *bl[2] = {"all vertices", "vertex #geq 10 tracks"};
    for (int a = 0; a < 2; ++a)
      for (int b = 0; b < 2; ++b)
      {
        gR[a][b]->SetMarkerStyle(ms[a][b]);
        gR[a][b]->SetMarkerColor(col[b]);
        gR[a][b]->SetLineColor(col[b]);
        gR[a][b]->SetMarkerSize(1.4);
        gR[a][b]->Draw("p same");
        lg->AddEntry(gR[a][b], Form("%s, %s", al[a], bl[b]), "p");
      }
    lg->Draw();
    tagline();
    save(c, "fig_trk_dcares");
  }

  // ntuple branch sizes table
  {
  std::map<std::string, std::string> brDesc = {
    {"run", "run number (EventHeader)"},
    {"event", "event sequence number (EventHeader)"},
    {"GL1_BCO", "GL1 beam-clock counter of the trigger"},
    {"GL1_scaled", "GL1 scaled trigger bit vector"},
    {"GL1_live", "GL1 live trigger bit vector"},
    {"GL1_bunch", "GL1 bunch number"},
    {"nTrk", "tracks in SvtxTrackMap"},
    {"nTrkSil", "tracks with $\\geq1$ MVTX or INTT cluster"},
    {"nTrkINTT", "tracks with $\\geq1$ INTT cluster"},
    {"nTrkMVTXonly", "MVTX but no INTT clusters"},
    {"nTrkTPConly", "no silicon cluster"},
    {"nTrkX0", "tracks with crossing 0 (incl.\\ TPC-only)"},
    {"nTrkSilX0", "silicon tracks with crossing 0"},
    {"nXingTrk", "distinct crossings among silicon tracks"},
    {"nVtx", "vertices in SvtxVertexMap"},
    {"nVtxX0", "vertices with crossing 0"},
    {"vx", "$x$ of the crossing-0 vertex with most tracks (cm, $-999$ if none)"},
    {"vy", "same, $y$ (cm)"},
    {"vz", "same, $z$ (cm)"},
    {"vntrk", "tracks in that vertex"},
    {"vchi2ndf", "its $\\chi^2$/ndf"},
    {"ievt", "entry number of the event in tree \\fcmd{event}"},
    {"id", "track key in SvtxTrackMap"},
    {"pt", "$p_T$ (GeV/$c$)"},
    {"eta", "pseudorapidity"},
    {"phi", "azimuth (rad)"},
    {"charge", "charge"},
    {"pcax", "track position $x$ at the PCA stored on the DST (cm)"},
    {"pcay", "same, $y$ (cm)"},
    {"pcaz", "same, $z$ (cm)"},
    {"chi2", "fit $\\chi^2$"},
    {"ndf", "fit ndf"},
    {"quality", "SvtxTrack quality ($\\chi^2$/ndf)"},
    {"crossing", "beam crossing relative to the GL1 BCO (0 for TPC-only)"},
    {"nmvtx", "MVTX clusters"},
    {"nintt", "INTT clusters"},
    {"ntpc", "TPC clusters"},
    {"ntpot", "TPOT clusters"},
    {"hasSiSeed", "track has a silicon seed (pointer or id)"},
    {"hasTpcSeed", "track has a TPC seed"},
    {"vtxMatch", "vertex association class 0--3 (Sec.~\\ref{sec:vtxassoc})"},
    {"vtxId", "key of the associated vertex ($-1$ none)"},
    {"tvx", "associated vertex $x$ (cm)"},
    {"tvy", "associated vertex $y$ (cm)"},
    {"tvz", "associated vertex $z$ (cm)"},
    {"tvntrk", "tracks in the associated vertex"},
    {"dcaxy", "signed transverse DCA to the associated vertex (cm)"},
    {"dcaxyErr", "its uncertainty (cm)"},
    {"dcaz", "longitudinal DCA to the associated vertex (cm)"},
    {"dcazErr", "its uncertainty (cm)"}};
    std::ofstream bt(note + "/tables/tab_branches.tex");
    bt << "% generated by scripts/make_qa.C\n";
    auto dump = [&](TTree *tr, long nev) {
      long tot = tr->GetZipBytes();
      for (auto *ob : *tr->GetListOfBranches())
      {
        auto *br = (TBranch *) ob;
        std::string n = br->GetName();
        std::string leaf = br->GetTitle();
        std::string ty = leaf.substr(leaf.find('/') + 1);
        bt << "\\fcmd{" << n << "} & " << ty << " & " << brDesc[n] << " & " << Form("%.1f", 100. * br->GetZipBytes() / tot) << "\\\\\n";
      }
      txt << tr->GetName() << " zip bytes " << tot << " (" << tot / 1024. / nev << " kB/event)\n";
      return tot;
    };
    bt << "\\midrule\n\\multicolumn{4}{l}{\\textbf{tree \\fcmd{event}} (one entry per event)}\\\\\n\\midrule\n";
    long ze = dump(te, nEv);
    bt << "\\midrule\n\\multicolumn{4}{l}{\\textbf{tree \\fcmd{track}} (one entry per track; also \\fcmd{run}, \\fcmd{event})}\\\\\n\\midrule\n";
    long zt = dump(tt, nEv);
    numf("evtTreekBPerEv", ze / 1024. / nEv, 2);
    numf("trkTreekBPerEv", zt / 1024. / nEv, 2);
  }
  // DST branch table
  {
    std::ofstream dt(note + "/tables/tab_dstsize.tex");
    dt << "% generated by scripts/make_qa.C\n\\begin{tabular}{lrr}\n\\toprule\nnode (branch of tree \\fcmd{T}) & DST (kB/ev) & trimmed DST (kB/ev)\\\\\n\\midrule\n";
    for (auto &[k, v] : dstBr)
    {
      std::string kk = std::regex_replace(std::regex_replace(k, std::regex("#"), "\\#"), std::regex("_"), "\\_");
      dt << "\\texttt{" << kk << "} & " << Form("%.2f", v / 1024. / nEv) << " & ";
      if (trimBr.count(k)) dt << Form("%.2f", trimBr[k] / 1024. / nEv);
      else dt << "dropped";
      dt << "\\\\\n";
    }
    dt << "\\bottomrule\n\\end{tabular}\n";
  }
  num("nChecks", std::to_string(nCheck));
  num("nFailed", std::to_string(nFail));
  txt << "FAILED CHECKS: " << nFail << "\n";
  tex.close();
  txt.close();
  // defaults file: \newcommand for every macro name (so main.tex compiles even before this script)
  {
    std::ifstream in(note + "/tables/qa_numbers.tex");
    std::ofstream out(note + "/tables/qa_defaults.tex");
    out << "% generated by scripts/make_qa.C: fallbacks for tables/qa_numbers.tex\n";
    std::string l;
    std::smatch m;
    while (std::getline(in, l))
      if (std::regex_search(l, m, std::regex("\\\\renewcommand\\{\\\\([A-Za-z]+)\\}"))) out << "\\newcommand{\\" << m[1] << "}{??}\n";
  }
  std::cout << "make_qa: failed checks " << nFail << std::endl;
}
