// fillV0MassCorr.C -- V0 mass vs pT histograms from the full KFP Ks / Lambda ntuples (covariance-free run 79509
// seg 0-30) with a per-daughter pT re-calibration from a kappa look-up table (MomentumScaleMap).
//
//   pT_corr = kappa(q, pT, eta, phi) * pT, eta and phi unchanged  ->  p_vec scaled by kappa (px, py, pz)
//   kappa: 4D table (q = +-1; pT 112 nodes; eta 13 nodes -1.5..1.5; phi 25 nodes -pi..pi), linear interpolation
//   in (ln pT, eta, phi), clamped at the table edges (the table itself is constant in ln pT beyond its clamp).
//
// Same selection as run3pp/macro/loosecov/fill_v0_mass.C (|PVz|<10, |eta_d|<1.1; QA: DIRA>0.995, Lxy>0.5,
// DCA12<0.2, |DCA_xy(d)|>0.02); the cuts use the stored KFP quantities (not recomputed after the correction).
// Daughter species from the KFP convention (Ks: pi pi; Lambda: track_2 = proton), PDG masses.
// mode: 0 = KFP fitted mass and pT (as fill_v0_mass.C), 1 = recomputed from daughters, no correction,
//       2 = recomputed from daughters with the kappa correction.
// Output histogram names as fill_v0_mass.C (hKs_qa, hLa_qa, hLaBar_qa, ...) so fitMomRes.C reads them directly.
//   root -l -b -q 'macros/fillV0MassCorr.C+("<Ks ntuple>","<Lambda ntuple>","<kappa csv>","<out.root>",mode)'
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TH2F.h>
#include <TProfile.h>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include <array>

namespace fc {
const double mPi = 0.13957039, mP = 0.93827209;

struct Kappa {
  std::vector<double> lpt, eta, phi;   // node coordinates (ln pT)
  std::vector<double> v;               // [iq][ipt][ieta][iphi]
  bool load(const char* fn)
  {
    std::ifstream in(fn);
    if (!in) { std::cerr << "cannot open " << fn << std::endl; return false; }
    std::string line;
    std::vector<std::array<double, 5>> rows;
    std::map<double, int> mpt, meta, mphi;
    bool header = false;
    while (std::getline(in, line)) {
      if (line.empty() || line[0] == '#') continue;
      if (!header) { header = true; continue; }   // q,pT,curv,eta,phi,kappa,...
      std::stringstream ss(line); std::string tok; std::vector<double> x;
      while (std::getline(ss, tok, ',')) x.push_back(std::stod(tok));
      if (x.size() < 6) continue;
      const double ph = std::round(x[4] * 1e5) / 1e5, et = std::round(x[3] * 1e5) / 1e5;
      rows.push_back({x[0], x[1], et, ph, x[5]});
      mpt[x[1]] = 0; meta[et] = 0; mphi[ph] = 0;
    }
    auto idx = [](std::map<double, int>& m, std::vector<double>& out, bool log) {
      int i = 0; for (auto& p : m) { p.second = i++; out.push_back(log ? std::log(p.first) : p.first); } };
    idx(mpt, lpt, true); idx(meta, eta, false); idx(mphi, phi, false);
    const size_t n = 2 * lpt.size() * eta.size() * phi.size();
    v.assign(n, NAN);
    for (auto& r : rows) {
      const int iq = r[0] > 0 ? 1 : 0;
      v[((iq * lpt.size() + mpt[r[1]]) * eta.size() + meta[r[2]]) * phi.size() + mphi[r[3]]] = r[4];
    }
    const long nan = std::count_if(v.begin(), v.end(), [](double x) { return std::isnan(x); });
    printf("kappa table: %zu rows, nodes pT %zu [%g,%g] eta %zu [%g,%g] phi %zu [%g,%g], missing %ld\n", rows.size(),
           lpt.size(), std::exp(lpt.front()), std::exp(lpt.back()), eta.size(), eta.front(), eta.back(), phi.size(),
           phi.front(), phi.back(), nan);
    return nan == 0;
  }
  static void locate(const std::vector<double>& g, double x, int& i, double& f)
  {
    if (x <= g.front()) { i = 0; f = 0; return; }
    if (x >= g.back()) { i = g.size() - 2; f = 1; return; }
    i = int(std::upper_bound(g.begin(), g.end(), x) - g.begin()) - 1;
    f = (x - g[i]) / (g[i + 1] - g[i]);
  }
  double operator()(int q, double pt, double et, double ph) const
  {
    while (ph < -M_PI) ph += 2 * M_PI;
    while (ph > M_PI) ph -= 2 * M_PI;
    int i, j, k; double fi, fj, fk;
    locate(lpt, std::log(pt), i, fi); locate(eta, et, j, fj); locate(phi, ph, k, fk);
    const int iq = q > 0 ? 1 : 0;
    const size_t ne = eta.size(), np = phi.size(), nt = lpt.size();
    double s = 0;
    for (int a = 0; a < 2; ++a)
      for (int b = 0; b < 2; ++b)
        for (int c = 0; c < 2; ++c) {
          const double w = (a ? fi : 1 - fi) * (b ? fj : 1 - fj) * (c ? fk : 1 - fk);
          if (w == 0) continue;
          s += w * v[((iq * nt + i + a) * ne + j + b) * np + k + c];
        }
    return s;
  }
};

long fill(TTree* t, const std::string& M, bool lambda, const Kappa& K, int mode, TH2F* hQa, TH2F* hLoose, TH2F* hQaAnti,
          TH2F* hLooseAnti, TProfile* pK1, TProfile* pK2, TH1D* hDm)
{
  float mass, dira, pt, lxy, d12, dxy1, dxy2, pvz, eta1, eta2, p1[3], p2[3];
  Char_t q1 = 0, q2 = 0;
  t->SetBranchStatus("*", 0);
  auto B = [&](const std::string& n, void* a) { t->SetBranchStatus(n.c_str(), 1); t->SetBranchAddress(n.c_str(), a); };
  B(M + "_mass", &mass); B(M + "_DIRA", &dira); B(M + "_pT", &pt); B(M + "_decayLength_xy", &lxy);
  B("track_1_track_2_DCA", &d12); B("track_1_PV_DCA_xy", &dxy1); B("track_2_PV_DCA_xy", &dxy2);
  B("primary_vertex_z", &pvz); B("track_1_pseudorapidity", &eta1); B("track_2_pseudorapidity", &eta2);
  B("track_1_px", &p1[0]); B("track_1_py", &p1[1]); B("track_1_pz", &p1[2]);
  B("track_2_px", &p2[0]); B("track_2_py", &p2[1]); B("track_2_pz", &p2[2]);
  B("track_1_charge", &q1); B("track_2_charge", &q2);
  const double m1 = lambda ? mPi : mPi, m2 = lambda ? mP : mPi;   // track_2 = proton for Lambda
  long n = 0;
  for (Long64_t i = 0; i < t->GetEntries(); ++i) {
    t->GetEntry(i);
    if (std::fabs(pvz) > 10 || std::fabs(eta1) > 1.1 || std::fabs(eta2) > 1.1) continue;
    double M0 = mass, PT = pt;
    if (mode >= 1) {
      double k1 = 1, k2 = 1;
      const double pt1 = std::hypot(p1[0], p1[1]), pt2 = std::hypot(p2[0], p2[1]);
      if (mode == 2) {
        k1 = K(q1, pt1, std::asinh(p1[2] / pt1), std::atan2(p1[1], p1[0]));
        k2 = K(q2, pt2, std::asinh(p2[2] / pt2), std::atan2(p2[1], p2[0]));
        pK1->Fill(pt1, k1); pK2->Fill(pt2, k2);
      }
      const double x1 = k1 * p1[0], y1 = k1 * p1[1], z1 = k1 * p1[2], x2 = k2 * p2[0], y2 = k2 * p2[1], z2 = k2 * p2[2];
      const double e1 = std::sqrt(x1 * x1 + y1 * y1 + z1 * z1 + m1 * m1), e2 = std::sqrt(x2 * x2 + y2 * y2 + z2 * z2 + m2 * m2);
      const double px = x1 + x2, py = y1 + y2, pz = z1 + z2, E = e1 + e2;
      M0 = std::sqrt(std::max(E * E - px * px - py * py - pz * pz, 0.));
      PT = std::hypot(px, py);
      hDm->Fill(M0 - mass);
    }
    // Lambda vs anti-Lambda from the proton (track_2) charge
    const bool anti = hQaAnti && q2 < 0;
    (anti ? hLooseAnti : hLoose)->Fill(PT, M0);
    if (dira < 0.995 || lxy < 0.5 || d12 > 0.2 || std::fabs(dxy1) < 0.02 || std::fabs(dxy2) < 0.02) continue;
    (anti ? hQaAnti : hQa)->Fill(PT, M0);
    ++n;
  }
  return n;
}
}  // namespace fc

void fillV0MassCorr(const char* ksFile, const char* laFile, const char* kappaCsv, const char* outFile, int mode = 2)
{
  fc::Kappa K;
  if (mode == 2 && !K.load(kappaCsv)) { std::cerr << "kappa table not usable" << std::endl; return; }
  if (mode == 2) {
    // spot checks against the csv nodes
    printf("kappa(q=-1,pT=1,eta=0,phi=0) = %.6f, kappa(q=+1,pT=1,eta=0,phi=0) = %.6f\n", K(-1, 1, 0, 0), K(1, 1, 0, 0));
  }
  TFile out(outFile, "RECREATE");
  TH2F* ksQa = new TH2F("hKs_qa", "K0S, QA cuts;p_{T} (GeV/c);M (GeV/c^{2})", 100, 0, 5, 200, 0.40, 0.60);
  TH2F* ksLo = new TH2F("hKs_loose", "K0S loose;p_{T} (GeV/c);M (GeV/c^{2})", 100, 0, 5, 200, 0.40, 0.60);
  TH2F* laQa = new TH2F("hLa_qa", "Lambda, QA cuts;p_{T};M", 100, 0, 5, 320, 1.08, 1.16);
  TH2F* laLo = new TH2F("hLa_loose", "Lambda, loose;p_{T};M", 100, 0, 5, 320, 1.08, 1.16);
  TH2F* lbQa = new TH2F("hLaBar_qa", "anti-Lambda, QA cuts;p_{T};M", 100, 0, 5, 320, 1.08, 1.16);
  TH2F* lbLo = new TH2F("hLaBar_loose", "anti-Lambda, loose;p_{T};M", 100, 0, 5, 320, 1.08, 1.16);
  TProfile* kK1 = new TProfile("pKappa_Ks_d1", "Ks track_1 (pi-);p_{T};kappa", 50, 0, 5);
  TProfile* kK2 = new TProfile("pKappa_Ks_d2", "Ks track_2 (pi+);p_{T};kappa", 50, 0, 5);
  TProfile* kL1 = new TProfile("pKappa_La_d1", "Lambda track_1 (pi);p_{T};kappa", 50, 0, 5);
  TProfile* kL2 = new TProfile("pKappa_La_d2", "Lambda track_2 (p);p_{T};kappa", 50, 0, 5);
  TH1D* dmK = new TH1D("hDm_Ks", "Ks: M(recomputed) - M(KFP);#DeltaM (GeV/c^{2})", 400, -0.02, 0.02);
  TH1D* dmL = new TH1D("hDm_La", "Lambda: M(recomputed) - M(KFP);#DeltaM (GeV/c^{2})", 400, -0.02, 0.02);
  TFile fk(ksFile), fl(laFile);
  TTree* tk = (TTree*)fk.Get("DecayTree");
  TTree* tl = (TTree*)fl.Get("DecayTree");
  out.cd();
  const long nk = fc::fill(tk, "K_S0", false, K, mode, ksQa, ksLo, nullptr, nullptr, kK1, kK2, dmK);
  const long nl = fc::fill(tl, "Lambda0", true, K, mode, laQa, laLo, lbQa, lbLo, kL1, kL2, dmL);
  out.cd();
  out.Write();
  out.Close();
  std::cout << "mode " << mode << ": Ks QA candidates " << nk << ", Lambda QA candidates " << nl << std::endl;
  std::cout << "All done" << std::endl;
}
