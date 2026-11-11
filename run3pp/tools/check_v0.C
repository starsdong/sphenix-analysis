// check_v0.C : QA of KFParticle Ks / Lambda ntuples with emphasis on the bunch-crossing choice.
//
// usage: root -b -q 'check_v0.C("outdir", "Ks")'   or  "Lambda"
//   reads  outdir/KFP_<ch>_*.root           (daughters same crossing + PV crossing match)
//          outdir/KFP_<ch>_noXmatch_*.root  (no crossing requirement at all)
//          outdir/KFP_EventInfo_*.root      (number of events)
//   writes outdir/qa/<ch>_*.png and prints a summary table.
//
// daughter type: I = has INTT clusters (crossing from INTT timing), M = MVTX but no INTT,
//                T = TPC-only (crossing not measurable, set to 0 by the tracking)
#include <TCanvas.h>
#include <TChain.h>
#include <TF1.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>

#include <cmath>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <tuple>
#include <vector>

namespace
{
  struct Cfg
  {
    std::string mother;
    double lo, hi;     // histogram range
    double m0;         // PDG mass
    double fitlo, fithi;
    int nbins;
  };

  char dtype(unsigned nintt, unsigned nmvtx) { return nintt > 0 ? 'I' : (nmvtx > 0 ? 'M' : 'T'); }
  std::string pairType(char a, char b)
  {
    std::string order = "IMT";
    if (order.find(a) > order.find(b)) std::swap(a, b);
    return std::string(1, a) + b;
  }

  struct FitRes
  {
    double n = 0, mean = 0, sigma = 0, emean = 0, esigma = 0, S = 0, B = 0;
    bool ok = false;
  };

  FitRes fitPeak(TH1F *h, const Cfg &c, int color)
  {
    FitRes r;
    r.n = h->GetEntries();
    if (h->Integral() < 30) return r;
    TF1 f(Form("f_%s", h->GetName()), "gaus(0)+pol2(3)", c.fitlo, c.fithi);
    double bin = h->GetBinWidth(1);
    double b0 = 0.5 * (h->GetBinContent(h->FindBin(c.fitlo + 2 * bin)) + h->GetBinContent(h->FindBin(c.fithi - 2 * bin)));
    double peak = std::max(1.0, h->GetBinContent(h->FindBin(c.m0)) - b0);
    double s0 = (c.m0 < 1 ? 0.008 : 0.003);
    f.SetParameters(peak, c.m0, s0, b0, 0, 0);
    f.SetParLimits(1, c.m0 - 3 * s0, c.m0 + 3 * s0);
    f.SetParLimits(2, 0.2 * s0, 5 * s0);
    f.SetParLimits(0, 0, 1e9);
    f.SetLineColor(color);
    f.SetLineWidth(2);
    int st = h->Fit(&f, "QRLS0");
    r.mean = f.GetParameter(1);
    r.sigma = std::fabs(f.GetParameter(2));
    r.emean = f.GetParError(1);
    r.esigma = f.GetParError(2);
    double lo = r.mean - 3 * r.sigma, hi = r.mean + 3 * r.sigma;
    TF1 g("g", "gaus", c.fitlo, c.fithi);
    g.SetParameters(f.GetParameter(0), f.GetParameter(1), f.GetParameter(2));
    TF1 p("p", "pol2", c.fitlo, c.fithi);
    p.SetParameters(f.GetParameter(3), f.GetParameter(4), f.GetParameter(5));
    r.S = g.Integral(lo, hi) / bin;
    r.B = p.Integral(lo, hi) / bin;
    r.ok = (st == 0);
    h->GetListOfFunctions()->Add(f.Clone());
    return r;
  }

  void printRes(const std::string &label, const FitRes &r)
  {
    printf("  %-34s N=%7.0f  mean=%8.4f+-%6.4f  sigma=%6.2f+-%5.2f MeV  S(3sig)=%7.1f  B=%7.1f  S/B=%5.2f %s\n",
           label.c_str(), r.n, r.mean, r.emean, 1e3 * r.sigma, 1e3 * r.esigma, r.S, r.B, r.B > 0 ? r.S / r.B : 0., r.ok ? "" : "(fit not converged/low stat)");
  }
}  // namespace

void check_v0(const std::string &dir = "output/test1", const std::string &ch = "Ks",
              double cutDIRA = 0.995, double cutDLxy = 0.5, double cutDCA = 0.2, double cutTrkDCAxy = 0.02)
{
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(1);
  Cfg c;
  if (ch == "Ks")
    c = {"K_S0", 0.40, 0.60, 0.497611, 0.42, 0.58, 100};
  else
    c = {"Lambda0", 1.05, 1.20, 1.115683, 1.08, 1.17, 150};
  const std::string M = c.mother;
  gSystem->mkdir((dir + "/qa").c_str(), true);

  TChain ev("EventInfo");
  ev.Add((dir + "/KFP_EventInfo_*.root").c_str());
  const double nEv = ev.GetEntries();
  printf("\n===== %s : %s   events analysed = %.0f   QA cuts: DIRA>%.3f, L_xy>%.2f cm, DCA(d1,d2)<%.2f cm, daughter DCA_xy(PV)>%.3f cm\n",
         ch.c_str(), dir.c_str(), nEv, cutDIRA, cutDLxy, cutDCA, cutTrkDCAxy);

  // ------------------------------------------------------------------ loop helper
  struct Cand
  {
    float mass, dira, dlxy, dl, dlerr, dca, d1xy, d2xy, x, y, z, pt;
    int x1, x2, xpv, q;
    char t1, t2;
    Long64_t cbco, gl1;
    int evt;
    float p1x, p2x;
  };
  auto readAll = [&](const std::string &pattern) {
    std::vector<Cand> v;
    TChain t("DecayTree");
    t.Add(pattern.c_str());
    TTreeReader r(&t);
    TTreeReaderValue<float> mass(r, (M + "_mass").c_str()), dira(r, (M + "_DIRA").c_str()),
        dlxy(r, (M + "_decayLength_xy").c_str()), dl(r, (M + "_decayLength").c_str()), dlerr(r, (M + "_decayLengthErr").c_str()),
        x(r, (M + "_x").c_str()), y(r, (M + "_y").c_str()), z(r, (M + "_z").c_str()), pt(r, (M + "_pT").c_str()),
        dca(r, "track_1_track_2_DCA"), d1xy(r, "track_1_PV_DCA_xy"), d2xy(r, "track_2_PV_DCA_xy"),
        p1x(r, "track_1_px"), p2x(r, "track_2_px");
    TTreeReaderValue<int> x1(r, "track_1_bunch_crossing"), x2(r, "track_2_bunch_crossing"), xpv(r, "primary_vertex_ID"),
        pdg1(r, "track_1_PDG_ID"), pdg2(r, "track_2_PDG_ID"), evt(r, "eventNumber");
    TTreeReaderValue<unsigned> i1(r, "track_1_INTT_nHits"), m1(r, "track_1_MVTX_nHits"), i2(r, "track_2_INTT_nHits"), m2(r, "track_2_MVTX_nHits");
    TTreeReaderValue<Long64_t> cbco(r, "Collision_BCO");
    TTreeReaderValue<Long64_t> gl1(r, "GL1_BCO");
    while (r.Next())
    {
      Cand k;
      k.mass = *mass; k.dira = *dira; k.dlxy = *dlxy; k.dl = *dl; k.dlerr = *dlerr; k.dca = *dca;
      k.d1xy = std::fabs(*d1xy); k.d2xy = std::fabs(*d2xy); k.x = *x; k.y = *y; k.z = *z; k.pt = *pt;
      k.x1 = *x1; k.x2 = *x2; k.xpv = *xpv;
      k.t1 = dtype(*i1, *m1); k.t2 = dtype(*i2, *m2);
      // particle / antiparticle tag (Lambda: sign of the proton)
      k.q = 0;
      if (std::abs(*pdg1) == 2212) k.q = (*pdg1 > 0 ? 1 : -1);
      if (std::abs(*pdg2) == 2212) k.q = (*pdg2 > 0 ? 1 : -1);
      k.cbco = *cbco; k.gl1 = *gl1; k.evt = *evt; k.p1x = *p1x; k.p2x = *p2x;
      v.push_back(k);
    }
    return v;
  };
  auto pass = [&](const Cand &k) {
    return k.dira > cutDIRA && k.dlxy > cutDLxy && k.dca < cutDCA && k.d1xy > cutTrkDCAxy && k.d2xy > cutTrkDCAxy;
  };

  // ================================================================== 1) default (crossing-matched) sample
  std::vector<Cand> mt = readAll(dir + "/KFP_" + ch + "_0*.root");
  printf("[matched] candidates in ntuple: %zu  (%.2f / event)\n", mt.size(), mt.size() / std::max(1.0, nEv));

  // duplicates: same physical candidate reconstructed in two overlapping streaming events
  //   -> identical Collision_BCO (= GL1 BCO + crossing) and identical daughter momenta
  std::map<std::tuple<Long64_t, int, int>, int> seen;
  int nDup = 0, nDupPass = 0;
  std::vector<bool> isDup(mt.size(), false);
  for (size_t i = 0; i < mt.size(); ++i)
  {
    auto key = std::make_tuple(mt[i].cbco, (int) std::lround(mt[i].mass * 1e5), (int) std::lround((mt[i].p1x + mt[i].p2x) * 1e4));
    if (seen.count(key))
    {
      isDup[i] = true;
      ++nDup;
      if (pass(mt[i])) ++nDupPass;
    }
    else
      seen[key] = i;
  }
  printf("[matched] duplicates across events (same Collision_BCO & kinematics): %d of %zu (%.1f%%); after QA cuts %d\n",
         nDup, mt.size(), 100. * nDup / std::max<size_t>(1, mt.size()), nDupPass);

  // consistency checks of the crossing logic
  int nBad = 0, nTT = 0;
  for (auto &k : mt)
  {
    if (k.x1 != k.x2 || k.x1 != k.xpv) ++nBad;
    if (k.t1 == 'T' && k.t2 == 'T') ++nTT;
  }
  printf("[matched] candidates violating X1==X2==X_PV: %d   ;  with two TPC-only daughters: %d\n", nBad, nTT);

  const std::vector<std::string> ptypes = {"II", "IM", "MM", "IT", "MT", "TT"};
  std::map<std::string, TH1F *> hm;
  TH1F *hAll = new TH1F(Form("h_%s_all", ch.c_str()), Form("%s matched, QA cuts, de-duplicated;m (GeV/c^{2});counts", ch.c_str()), c.nbins, c.lo, c.hi);
  TH1F *hAllDup = new TH1F(Form("h_%s_alldup", ch.c_str()), "", c.nbins, c.lo, c.hi);
  TH1F *hPart = new TH1F(Form("h_%s_part", ch.c_str()), "", c.nbins, c.lo, c.hi);
  TH1F *hAnti = new TH1F(Form("h_%s_anti", ch.c_str()), "", c.nbins, c.lo, c.hi);
  for (auto &p : ptypes) hm[p] = new TH1F(Form("h_%s_%s", ch.c_str(), p.c_str()), Form("%s;m (GeV/c^{2});counts", p.c_str()), c.nbins, c.lo, c.hi);
  TH1F *hX = new TH1F(Form("hX_%s", ch.c_str()), Form("%s matched: daughter crossing;crossing (relative to GL1 BCO);candidates", ch.c_str()), 200, -150, 450);
  TH1F *hR = new TH1F(Form("hR_%s", ch.c_str()), Form("%s matched, |m-m_{0}|<2#sigma: decay radius;r_{xy} (cm);candidates", ch.c_str()), 60, 0, 30);
  std::map<std::string, TH1F *> hRt;
  for (auto &p : ptypes) hRt[p] = new TH1F(Form("hR_%s_%s", ch.c_str(), p.c_str()), "", 60, 0, 30);

  for (size_t i = 0; i < mt.size(); ++i)
  {
    auto &k = mt[i];
    if (!pass(k)) continue;
    hAllDup->Fill(k.mass);
    if (isDup[i]) continue;
    hAll->Fill(k.mass);
    hm[pairType(k.t1, k.t2)]->Fill(k.mass);
    hX->Fill(k.x1);
    if (ch == "Lambda") (k.q > 0 ? hPart : hAnti)->Fill(k.mass);
  }
  printf("[matched] after QA cuts: %.0f candidates incl. duplicates, %.0f de-duplicated\n", hAllDup->GetEntries(), hAll->GetEntries());
  FitRes rAll = fitPeak(hAll, c, kRed);
  printRes("ALL (de-dup)", rAll);
  if (ch == "Lambda")
  {
    printRes("  Lambda    (p+)", fitPeak(hPart, c, kRed));
    printRes("  antiLambda(p-)", fitPeak(hAnti, c, kRed));
  }
  for (auto &p : ptypes) printRes("  daughters " + p, fitPeak(hm[p], c, kRed));
  // decay radius of peak candidates per pair type
  for (size_t i = 0; i < mt.size(); ++i)
  {
    auto &k = mt[i];
    if (!pass(k) || isDup[i] || !rAll.ok) continue;
    if (std::fabs(k.mass - rAll.mean) < 2 * rAll.sigma)
    {
      double r = std::hypot(k.x, k.y);
      hR->Fill(r);
      hRt[pairType(k.t1, k.t2)]->Fill(r);
    }
  }

  // ================================================================== 2) no crossing requirement
  std::vector<Cand> nx = readAll(dir + "/KFP_" + ch + "_noXmatch_*.root");
  printf("[noXmatch] candidates in ntuple: %zu  (%.2f / event)\n", nx.size(), nx.size() / std::max(1.0, nEv));
  // classes
  //  A  X1==X2==Xpv              (what the default selection keeps)
  //  B  X1==X2 != Xpv            (pair consistent, PV from another crossing)
  //  C  X1!=X2, both silicon      (pair from different crossings: random)
  //  D  one TPC-only (X=0) + one silicon with X!=0   (crossing of TPC-only track unknown)
  //  E  one TPC-only + silicon with X==0
  //  F  both TPC-only
  auto cls = [](const Cand &k) -> std::string {
    bool T1 = k.t1 == 'T', T2 = k.t2 == 'T';
    if (T1 && T2) return "F: both TPC-only";
    if (T1 || T2)
    {
      int xs = T1 ? k.x2 : k.x1;
      return xs == 0 ? "E: TPC-only + Si(X=0)" : "D: TPC-only + Si(X!=0)";
    }
    if (k.x1 == k.x2) return k.x1 == k.xpv ? "A: X1=X2=Xpv" : "B: X1=X2!=Xpv";
    return "C: X1!=X2 (both Si)";
  };
  const std::vector<std::string> classes = {"A: X1=X2=Xpv", "B: X1=X2!=Xpv", "C: X1!=X2 (both Si)", "D: TPC-only + Si(X!=0)", "E: TPC-only + Si(X=0)", "F: both TPC-only"};
  std::map<std::string, TH1F *> hc;
  int ic = 0;
  for (auto &cl : classes) hc[cl] = new TH1F(Form("hc_%s_%d", ch.c_str(), ic++), (cl + ";m (GeV/c^{2});counts").c_str(), c.nbins, c.lo, c.hi);
  TH1F *hdX = new TH1F(Form("hdX_%s", ch.c_str()), Form("%s noXmatch, both Si, QA cuts;X_{1}-X_{2};candidates", ch.c_str()), 201, -100.5, 100.5);
  TH1F *hdXpv = new TH1F(Form("hdXpv_%s", ch.c_str()), Form("%s noXmatch, X1=X2 (Si), QA cuts;X_{daughters}-X_{PV};candidates", ch.c_str()), 201, -100.5, 100.5);
  TH1F *hDpeak = new TH1F(Form("hDpeak_%s", ch.c_str()), Form("%s noXmatch class D in peak;X_{Si daughter};candidates", ch.c_str()), 120, -150, 450);
  std::set<std::tuple<int, int, int>> nxSeen;
  for (auto &k : nx)
  {
    if (!pass(k)) continue;
    std::string cl = cls(k);
    hc[cl]->Fill(k.mass);
    if (k.t1 != 'T' && k.t2 != 'T')
    {
      hdX->Fill(k.x1 - k.x2);
      if (k.x1 == k.x2) hdXpv->Fill(k.x1 - k.xpv);
    }
  }
  std::map<std::string, FitRes> rc;
  for (auto &cl : classes)
  {
    rc[cl] = fitPeak(hc[cl], c, kBlue);
    printRes("  " + cl, rc[cl]);
  }
  for (auto &k : nx)
  {
    if (!pass(k) || cls(k)[0] != 'D' || !rAll.ok) continue;
    if (std::fabs(k.mass - rAll.mean) < 2 * rAll.sigma) hDpeak->Fill(k.t1 == 'T' ? k.x2 : k.x1);
  }

  // ================================================================== 3) fit-free sideband subtraction
  // peak window |m-m_ref| < 3 sigma_ref, sidebands 6-9 sigma on both sides (same total width)
  double mref = rAll.ok ? rAll.mean : c.m0;
  double sref = (rAll.ok && rAll.sigma > 0.5e-3) ? rAll.sigma : (ch == "Ks" ? 0.005 : 0.0015);
  auto region = [&](float m) { double d = std::fabs(m - mref) / sref; return d < 3 ? 1 : (d > 6 && d < 9 ? -1 : 0); };
  auto sbs = [&](const std::vector<Cand> &v, const std::vector<bool> *skip, auto pred, double &err) {
    double np = 0, ns = 0;
    for (size_t i = 0; i < v.size(); ++i)
    {
      if (skip && (*skip)[i]) continue;
      const auto &k = v[i];
      if (!pass(k) || !pred(k)) continue;
      int rg = region(k.mass);
      if (rg == 1) ++np;
      if (rg == -1) ++ns;
    }
    err = std::sqrt(np + ns);
    return np - ns;
  };
  printf("--- sideband-subtracted signal (peak |m-%.4f|<3x%.1f MeV, sidebands 6-9 sigma), QA cuts\n", mref, 1e3 * sref);
  double e;
  printf("  matched, all (de-dup)              : %7.1f +- %5.1f\n", sbs(mt, &isDup, [](const Cand &) { return true; }, e), e);
  for (auto &p : ptypes)
  {
    double s = sbs(mt, &isDup, [&](const Cand &k) { return pairType(k.t1, k.t2) == p; }, e);
    printf("  matched, daughters %s              : %7.1f +- %5.1f\n", p.c_str(), s, e);
  }
  for (auto &cl : classes)
  {
    double s = sbs(nx, nullptr, [&](const Cand &k) { return cls(k) == cl; }, e);
    printf("  noXmatch %-26s: %7.1f +- %5.1f\n", cl.c_str(), s, e);
  }
  // sideband-subtracted Delta-crossing between the two silicon daughters (true V0s should have dX = 0)
  TH1F *hdXsig = new TH1F(Form("hdXsig_%s", ch.c_str()), Form("%s noXmatch both Si, peak - sidebands;X_{1}-X_{2};signal candidates", ch.c_str()), 81, -40.5, 40.5);
  TH1F *hdXpvsig = new TH1F(Form("hdXpvsig_%s", ch.c_str()), Form("%s noXmatch X1=X2 (Si), peak - sidebands;X_{daughters}-X_{PV};signal candidates", ch.c_str()), 81, -40.5, 40.5);
  for (auto &k : nx)
  {
    if (!pass(k) || k.t1 == 'T' || k.t2 == 'T') continue;
    int rg = region(k.mass);
    if (rg == 0) continue;
    hdXsig->Fill(std::max(-40, std::min(40, k.x1 - k.x2)), rg);
    if (k.x1 == k.x2) hdXpvsig->Fill(std::max(-40, std::min(40, k.x1 - k.xpv)), rg);
  }
  printf("  noXmatch both-Si signal: dX=0 %.1f, dX!=0 %.1f  (fraction of true V0s with daughters in different crossings: %.2f)\n",
         hdXsig->GetBinContent(hdXsig->FindBin(0)), hdXsig->Integral() - hdXsig->GetBinContent(hdXsig->FindBin(0)),
         (hdXsig->Integral() - hdXsig->GetBinContent(hdXsig->FindBin(0))) / std::max(1e-9, hdXsig->Integral()));
  // which daughter types appear among dX!=0 signal
  {
    std::map<std::string, double> m;
    for (auto &k : nx)
    {
      if (!pass(k) || k.t1 == 'T' || k.t2 == 'T' || k.x1 == k.x2) continue;
      int rg = region(k.mass);
      if (rg) m[pairType(k.t1, k.t2)] += rg;
    }
    printf("  noXmatch both-Si, dX!=0, signal by daughter type:");
    for (auto &[p, n] : m) printf("  %s %.1f", p.c_str(), n);
    printf("\n");
  }

  // ================================================================== plots
  const std::string pre = dir + "/qa/" + ch;
  TCanvas cv("cv", "", 1500, 1000);
  cv.Divide(3, 2);
  for (size_t i = 0; i < ptypes.size(); ++i)
  {
    cv.cd(i + 1);
    hm[ptypes[i]]->SetTitle(Form("%s matched, daughters %s (I=INTT, M=MVTX-only, T=TPC-only)", ch.c_str(), ptypes[i].c_str()));
    hm[ptypes[i]]->SetMarkerStyle(20);
    hm[ptypes[i]]->SetMarkerSize(0.6);
    hm[ptypes[i]]->Draw("E");
    if (auto *f = hm[ptypes[i]]->GetFunction(Form("f_%s", hm[ptypes[i]]->GetName()))) f->Draw("same");
  }
  cv.SaveAs((pre + "_matched_by_daughter_type.png").c_str());

  TCanvas cv2("cv2", "", 1500, 1000);
  cv2.Divide(3, 2);
  for (size_t i = 0; i < classes.size(); ++i)
  {
    cv2.cd(i + 1);
    hc[classes[i]]->SetTitle((ch + " no crossing req., " + classes[i]).c_str());
    hc[classes[i]]->SetMarkerStyle(20);
    hc[classes[i]]->SetMarkerSize(0.6);
    hc[classes[i]]->Draw("E");
    if (auto *f = hc[classes[i]]->GetFunction(Form("f_%s", hc[classes[i]]->GetName()))) f->Draw("same");
  }
  cv2.SaveAs((pre + "_noXmatch_by_crossing_class.png").c_str());

  TCanvas cv3("cv3", "", 1500, 1000);
  cv3.Divide(3, 2);
  cv3.cd(1);
  hAll->SetMarkerStyle(20);
  hAll->SetMarkerSize(0.6);
  hAll->Draw("E");
  if (auto *f = hAll->GetFunction(Form("f_%s", hAll->GetName()))) f->Draw("same");
  TLatex lt;
  lt.SetNDC();
  lt.SetTextSize(0.045);
  lt.DrawLatex(0.15, 0.85, Form("mean %.4f, #sigma %.1f MeV", rAll.mean, 1e3 * rAll.sigma));
  lt.DrawLatex(0.15, 0.79, Form("S=%.0f B=%.0f (3#sigma), N_{ev}=%.0f", rAll.S, rAll.B, nEv));
  cv3.cd(2);
  hX->Draw();
  cv3.cd(3);
  hR->SetLineColor(kBlack);
  hR->Draw();
  int col[6] = {kRed, kBlue, kGreen + 2, kMagenta, kOrange + 1, kCyan + 2};
  TLegend lg(0.6, 0.5, 0.88, 0.88);
  for (size_t i = 0; i < ptypes.size(); ++i)
  {
    hRt[ptypes[i]]->SetLineColor(col[i]);
    hRt[ptypes[i]]->Draw("same");
    lg.AddEntry(hRt[ptypes[i]], ptypes[i].c_str(), "l");
  }
  lg.Draw();
  cv3.cd(4);
  hdXsig->SetMarkerStyle(20);
  hdXsig->Draw("E");
  cv3.cd(5);
  hdXpvsig->SetMarkerStyle(20);
  hdXpvsig->Draw("E");
  cv3.cd(6);
  gPad->SetLogy();
  hdX->Draw();
  cv3.SaveAs((pre + "_crossing_overview.png").c_str());
}
