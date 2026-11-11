// fill_v0_mass_kappa.C -- reduce KFParticle Ks / Lambda ntuples (loosecov production) to invariant-mass vs pT
// histograms, in three versions, for the DDSim momentum calibration with the kappa per-track pT re-calibration.
//   root -l -b -q 'fill_v0_mass_kappa.C("<kfp dir>","<out prefix>","<kappa csv>")'
// Writes <out prefix>_kfp.root   : KFParticle fitted mass and pT (identical to fill_v0_mass.C histograms)
//        <out prefix>_raw.root   : mass and V0 pT recomputed from the stored daughter momenta (PDG masses), no correction
//        <out prefix>_kappa.root : same, daughter momenta scaled by kappa(q, pT, eta, phi) (eta, phi unchanged)
// Histogram names as fill_v0_mass.C (hKs_qa, hLa_qa, hLaBar_qa, *_loose) + nEventsProcessed, so fitMomRes.C reads any of them.
// Cuts (stored KFP values, not recomputed): all |PV z| < 10 cm, |eta_d| < 1.1; "qa" + DIRA > 0.995, L_xy > 0.5 cm,
// DCA12 < 0.2 cm, |DCA_xy(d)| > 0.02 cm. Daughters: Ks track_1 = pi-, track_2 = pi+; Lambda track_2 = (anti)proton,
// Lambda vs anti-Lambda from the proton charge (all three versions).
// kappa table: csv q,pT,curv,eta,phi,kappa,...; linear interpolation in (ln pT, eta, phi), clamped at the table edges.
// Same code as DataDrivenSim/macros/fillV0MassCorr.C (validated there on run 79509 seg 0-30).
#include <TFile.h>
#include <TChain.h>
#include <TH2F.h>
#include <TParameter.h>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <array>
#include <algorithm>

namespace fk {
const double mPi = 0.13957039, mP = 0.93827209;

struct Kappa {
  std::vector<double> lpt, eta, phi, v;
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
      if (!header) { header = true; continue; }
      std::stringstream ss(line); std::string tok; std::vector<double> x;
      while (std::getline(ss, tok, ',')) x.push_back(std::stod(tok));
      if (x.size() < 6) continue;
      const double ph = std::round(x[4] * 1e5) / 1e5, et = std::round(x[3] * 1e5) / 1e5;
      rows.push_back({x[0], x[1], et, ph, x[5]});
      mpt[x[1]] = 0; meta[et] = 0; mphi[ph] = 0;
    }
    auto idx = [](std::map<double, int>& m, std::vector<double>& out, bool lg) {
      int i = 0; for (auto& p : m) { p.second = i++; out.push_back(lg ? std::log(p.first) : p.first); } };
    idx(mpt, lpt, true); idx(meta, eta, false); idx(mphi, phi, false);
    v.assign(2 * lpt.size() * eta.size() * phi.size(), NAN);
    for (auto& r : rows) {
      const int iq = r[0] > 0 ? 1 : 0;
      v[((iq * lpt.size() + mpt[r[1]]) * eta.size() + meta[r[2]]) * phi.size() + mphi[r[3]]] = r[4];
    }
    const long nan = std::count_if(v.begin(), v.end(), [](double x) { return std::isnan(x); });
    printf("kappa table: %zu rows, nodes pT %zu eta %zu phi %zu, missing %ld\n", rows.size(), lpt.size(), eta.size(), phi.size(), nan);
    return nan == 0 && !rows.empty();
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
    const size_t ne = eta.size(), np = phi.size(), nt = lpt.size();
    const int iq = q > 0 ? 1 : 0;
    double s = 0;
    for (int a = 0; a < 2; ++a)
      for (int b = 0; b < 2; ++b)
        for (int c = 0; c < 2; ++c) {
          const double w = (a ? fi : 1 - fi) * (b ? fj : 1 - fj) * (c ? fk : 1 - fk);
          if (w > 0) s += w * v[((iq * nt + i + a) * ne + j + b) * np + k + c];
        }
    return s;
  }
};

struct Set {   // one version (kfp / raw / kappa)
  TH2F *ksQa, *ksLo, *laQa, *laLo, *lbQa, *lbLo;
  explicit Set(const char* t)
  {
    ksQa = new TH2F(Form("hKs_qa_%s", t), "K0S, QA cuts;p_{T} (GeV/c);M (GeV/c^{2})", 100, 0, 5, 200, 0.40, 0.60);
    ksLo = new TH2F(Form("hKs_loose_%s", t), "K0S, |PVz|<10 |eta|<1.1 only;p_{T} (GeV/c);M (GeV/c^{2})", 100, 0, 5, 200, 0.40, 0.60);
    laQa = new TH2F(Form("hLa_qa_%s", t), "Lambda, QA cuts;p_{T};M", 100, 0, 5, 320, 1.08, 1.16);
    laLo = new TH2F(Form("hLa_loose_%s", t), "Lambda, loose;p_{T};M", 100, 0, 5, 320, 1.08, 1.16);
    lbQa = new TH2F(Form("hLaBar_qa_%s", t), "anti-Lambda, QA cuts;p_{T};M", 100, 0, 5, 320, 1.08, 1.16);
    lbLo = new TH2F(Form("hLaBar_loose_%s", t), "anti-Lambda, loose;p_{T};M", 100, 0, 5, 320, 1.08, 1.16);
    for (TH2F* h : {ksQa, ksLo, laQa, laLo, lbQa, lbLo}) h->SetDirectory(nullptr);
  }
  void write(const char* fn, Long64_t nev)
  {
    TFile o(fn, "RECREATE");
    ksQa->Write("hKs_qa"); ksLo->Write("hKs_loose"); laQa->Write("hLa_qa"); laLo->Write("hLa_loose");
    lbQa->Write("hLaBar_qa"); lbLo->Write("hLaBar_loose");
    TParameter<Long64_t>("nEventsProcessed", nev).Write();
    o.Close();
  }
};

long fill(TChain& t, const std::string& M, bool lambda, const Kappa& K, Set* S[3], double sumK[4])
{
  float mass, dira, pt, lxy, d12, dxy1, dxy2, pvz, eta1, eta2, p1[3], p2[3];
  Char_t q1 = 0, q2 = 0;
  t.SetBranchStatus("*", 0);
  auto B = [&](const std::string& n, void* a) { t.SetBranchStatus(n.c_str(), 1); t.SetBranchAddress(n.c_str(), a); };
  B(M + "_mass", &mass); B(M + "_DIRA", &dira); B(M + "_pT", &pt); B(M + "_decayLength_xy", &lxy);
  B("track_1_track_2_DCA", &d12); B("track_1_PV_DCA_xy", &dxy1); B("track_2_PV_DCA_xy", &dxy2);
  B("primary_vertex_z", &pvz); B("track_1_pseudorapidity", &eta1); B("track_2_pseudorapidity", &eta2);
  B("track_1_px", &p1[0]); B("track_1_py", &p1[1]); B("track_1_pz", &p1[2]);
  B("track_2_px", &p2[0]); B("track_2_py", &p2[1]); B("track_2_pz", &p2[2]);
  B("track_1_charge", &q1); B("track_2_charge", &q2);
  const double m1 = mPi, m2 = lambda ? mP : mPi;
  long n = 0;
  for (Long64_t i = 0; i < t.GetEntries(); ++i) {
    t.GetEntry(i);
    if (std::fabs(pvz) > 10 || std::fabs(eta1) > 1.1 || std::fabs(eta2) > 1.1) continue;
    const bool qa = !(dira < 0.995 || lxy < 0.5 || d12 > 0.2 || std::fabs(dxy1) < 0.02 || std::fabs(dxy2) < 0.02);
    const double pt1 = std::hypot(p1[0], p1[1]), pt2 = std::hypot(p2[0], p2[1]);
    const double k1 = K(q1, pt1, std::asinh(p1[2] / pt1), std::atan2(p1[1], p1[0]));
    const double k2 = K(q2, pt2, std::asinh(p2[2] / pt2), std::atan2(p2[1], p2[0]));
    if (qa) { sumK[0] += k1; sumK[1] += k2; sumK[2] += 1; }
    const bool anti = lambda && q2 < 0;
    for (int v = 0; v < 3; ++v) {
      double M0 = mass, PT = pt;
      if (v > 0) {
        const double f1 = v == 2 ? k1 : 1, f2 = v == 2 ? k2 : 1;
        const double x1 = f1 * p1[0], y1 = f1 * p1[1], z1 = f1 * p1[2], x2 = f2 * p2[0], y2 = f2 * p2[1], z2 = f2 * p2[2];
        const double e1 = std::sqrt(x1 * x1 + y1 * y1 + z1 * z1 + m1 * m1), e2 = std::sqrt(x2 * x2 + y2 * y2 + z2 * z2 + m2 * m2);
        const double px = x1 + x2, py = y1 + y2, pz = z1 + z2, E = e1 + e2;
        M0 = std::sqrt(std::max(E * E - px * px - py * py - pz * pz, 0.));
        PT = std::hypot(px, py);
      }
      if (!lambda) { S[v]->ksLo->Fill(PT, M0); if (qa) S[v]->ksQa->Fill(PT, M0); }
      else { (anti ? S[v]->lbLo : S[v]->laLo)->Fill(PT, M0); if (qa) (anti ? S[v]->lbQa : S[v]->laQa)->Fill(PT, M0); }
    }
    if (qa) ++n;
  }
  return n;
}
}  // namespace fk

void fill_v0_mass_kappa(const char* kfpDir, const char* outPrefix, const char* kappaCsv)
{
  fk::Kappa K;
  if (!K.load(kappaCsv)) { std::cerr << "kappa table not usable" << std::endl; return; }
  printf("kappa(q=-1,pT=1,eta=0,phi=0) = %.6f, kappa(q=+1,pT=1,eta=0,phi=0) = %.6f\n", K(-1, 1, 0, 0), K(1, 1, 0, 0));
  TChain tk("DecayTree"), tl("DecayTree"), te("EventInfo");
  tk.Add(Form("%s/KFP_Ks_*.root", kfpDir));
  tl.Add(Form("%s/KFP_Lambda_*.root", kfpDir));
  te.Add(Form("%s/KFP_EventInfo_*.root", kfpDir));
  fk::Set* S[3] = {new fk::Set("kfp"), new fk::Set("raw"), new fk::Set("kappa")};
  double sk[4] = {0, 0, 0, 0}, sl[4] = {0, 0, 0, 0};
  const long nk = tk.GetNtrees() ? fk::fill(tk, "K_S0", false, K, S, sk) : 0;
  const long nl = tl.GetNtrees() ? fk::fill(tl, "Lambda0", true, K, S, sl) : 0;
  const Long64_t nev = te.GetEntries();
  const char* tag[3] = {"kfp", "raw", "kappa"};
  for (int v = 0; v < 3; ++v) S[v]->write(Form("%s_%s.root", outPrefix, tag[v]), nev);
  std::cout << "Ks QA candidates " << nk << ", Lambda QA candidates " << nl << ", events " << nev << std::endl;
  printf("mean kappa (QA): Ks d1 %.4f d2 %.4f ; Lambda d1 %.4f d2 %.4f\n", sk[2] ? sk[0] / sk[2] : 0, sk[2] ? sk[1] / sk[2] : 0,
         sl[2] ? sl[0] / sl[2] : 0, sl[2] ? sl[1] / sl[2] : 0);
  std::cout << "All done" << std::endl;
}
