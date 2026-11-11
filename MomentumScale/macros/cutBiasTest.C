// cutBiasTest.C -- does a tighter topological cut select better-measured candidates in data in a way DDSim cannot follow?
// DDSim applies the selection BEFORE the momentum smearing (directions are not smeared), so any narrowing of the data peak
// caused by the cut through resolution correlations is absent in the templates and would bias (a,b) low.
// Test: Gaussian width of the mass peak with the base QA cuts and with the tight cuts, in data (KFP ntuples, run 79509
// seg 0-30, KFP mass) and in DDSim (templates smeared at the main point a=0.0375, b=0.014, s=0.991, weights wLevy).
// Compare the ratio sigma(tight)/sigma(base) between data and DDSim per V0-pT bin.
//   root -l -b -q 'macros/cutBiasTest.C+(0.9995, 0.10, 0.03)'      args: DIRA min, |DCAxy(pi)| min, |DCAxy(p)| min (Lambda)
// Ks tight = DIRA min and |DCAxy(pi)| min for both pions. Prints BIAS rows.
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TF1.h>
#include <TRandom3.h>
#include <cmath>
#include <cstdio>
#include <vector>

namespace cb {
const double mPi = 0.13957, mP = 0.938272;
struct G { double s, es, m, em; };
G gfit(TH1D* h, bool ks)
{
  const double m0 = ks ? 0.4976 : 1.1157, lo = ks ? 0.46 : 1.103, hi = ks ? 0.54 : 1.130;
  TF1 f("f", ks ? "gaus(0)+pol1(3)" : "gaus(0)+pol2(3)", lo, hi);
  // initial background: straight line through the side-band averages
  auto avg = [&](double x1, double x2) { double t = 0; int n = 0; for (int i = h->FindBin(x1); i <= h->FindBin(x2); ++i) { t += h->GetBinContent(i); ++n; } return t / std::max(n, 1); };
  const double w = ks ? 0.012 : 0.004;
  const double yl = avg(lo, lo + w), yr = avg(hi - w, hi), sl = (yr - yl) / (hi - lo - w);
  const double c0 = yl - sl * (lo + 0.5 * w);
  f.SetParameters(std::max(h->GetBinContent(h->FindBin(m0)) - (c0 + sl * m0), 1.), m0, ks ? 0.005 : 0.0018, c0, sl, 0);
  if (!ks) { f.SetParameter(3, c0); f.SetParameter(4, sl); f.SetParameter(5, 0); }
  f.SetParLimits(1, m0 - 0.006, m0 + 0.006); f.SetParLimits(2, 0.0003, ks ? 0.03 : 0.012);
  h->Fit(&f, "QNRL");
  return {1e3 * std::fabs(f.GetParameter(2)), 1e3 * f.GetParError(2), 1e3 * f.GetParameter(1), 1e3 * f.GetParError(1)};
}
double smear(double p3[3], double m, double a, double b, double s, double g, double out[3])
{
  const double pt = std::hypot(p3[0], p3[1]), p = std::sqrt(pt * pt + p3[2] * p3[2]), beta = p / std::sqrt(p * p + m * m);
  const double sig = std::sqrt(std::pow(a * pt, 2) + std::pow(b / beta, 2));
  const double f = s * (1 + sig * g);
  for (int k = 0; k < 3; ++k) out[k] = f * p3[k];
  return 0;
}
double minv(double* a, double ma, double* b, double mb, double& pt)
{
  const double ea = std::sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2] + ma * ma), eb = std::sqrt(b[0] * b[0] + b[1] * b[1] + b[2] * b[2] + mb * mb);
  const double px = a[0] + b[0], py = a[1] + b[1], pz = a[2] + b[2];
  pt = std::hypot(px, py);
  return std::sqrt(std::max(0., (ea + eb) * (ea + eb) - px * px - py * py - pz * pz));
}
}  // namespace cb

void cutBiasTest(double diraMin = 0.9995, double piDca = 0.10, double pDca = 0.03)
{
  using namespace cb;
  for (int c = 0; c < 2; ++c) {
    const bool ks = c == 0;
    const std::vector<double> pb = ks ? std::vector<double>{0.2, 0.7, 1.2, 2.2, 4.0} : std::vector<double>{0.5, 1.2, 4.0};
    const int nb = pb.size() - 1;
    const int nbin = ks ? 60 : 120; const double hlo = ks ? 0.44 : 1.09, hhi = ks ? 0.56 : 1.15;
    std::vector<TH1D*> hd[2], hs[2];
    for (int t = 0; t < 2; ++t)
      for (int j = 0; j <= nb; ++j) {
        hd[t].push_back(new TH1D(Form("d%d%d%d", c, t, j), "", nbin, hlo, hhi));
        hs[t].push_back(new TH1D(Form("s%d%d%d", c, t, j), "", ks ? 240 : 240, hlo, hhi));
        hs[t].back()->Sumw2();
      }
    auto bin = [&](double pt) { if (pt < pb[0] || pt >= pb[nb]) return -1; int k = 0; while (pt >= pb[k + 1]) ++k; return k; };
    // ---- data
    {
      TFile f(ks ? "data/loosecov/KFP_Ks_run79509_seg000-030_loosecov.root" : "data/loosecov/KFP_Lambda_run79509_seg000-030_loosecov.root");
      TTree* t = (TTree*)f.Get("DecayTree");
      const char* M = ks ? "K_S0" : "Lambda0";
      float m, pt, dira, lxy, d12, pvz, e1, e2, x1, x2;
      t->SetBranchStatus("*", 0);
      std::vector<std::pair<std::string, float*>> br = {{Form("%s_mass", M), &m}, {Form("%s_pT", M), &pt}, {Form("%s_DIRA", M), &dira},
        {Form("%s_decayLength_xy", M), &lxy}, {"track_1_track_2_DCA", &d12}, {"primary_vertex_z", &pvz}, {"track_1_pseudorapidity", &e1},
        {"track_2_pseudorapidity", &e2}, {"track_1_PV_DCA_xy", &x1}, {"track_2_PV_DCA_xy", &x2}};
      for (auto& b : br) { t->SetBranchStatus(b.first.c_str(), 1); t->SetBranchAddress(b.first.c_str(), b.second); }
      for (Long64_t i = 0; i < t->GetEntries(); ++i) {
        t->GetEntry(i);
        if (std::fabs(pvz) >= 10 || std::fabs(e1) >= 1.1 || std::fabs(e2) >= 1.1) continue;
        if (!(dira > 0.995 && lxy > 0.5 && d12 < 0.2 && std::fabs(x1) > 0.02 && std::fabs(x2) > 0.02)) continue;
        const int k = bin(pt); if (k < 0) continue;
        // track_1 = pi- (Ks) / pion (Lambda); track_2 = pi+ (Ks) / proton (Lambda)
        const bool tight = dira > diraMin && std::fabs(x1) > piDca && std::fabs(x2) > (ks ? piDca : pDca);
        for (int tt = 0; tt < (tight ? 2 : 1); ++tt) { hd[tt][k]->Fill(m); hd[tt][nb]->Fill(m); }
      }
    }
    // ---- DDSim
    {
      TFile f(ks ? "output/momfit/Ks_nosmear_v2.root" : "output/momfit/Lambda_nosmear_v2.root");
      TTree* t = (TTree*)f.Get("ddsim");
      const char* d1 = ks ? "PiM" : "Pi"; const char* d2 = ks ? "PiP" : "P";
      float p1[3], p2[3], w, cosT, dx1, dx2; bool pass;
      const char* cc[3] = {"x", "y", "z"};
      t->SetBranchStatus("*", 0);
      for (int k = 0; k < 3; ++k) {
        t->SetBranchStatus(Form("rc%sSvP%s", d1, cc[k]), 1); t->SetBranchAddress(Form("rc%sSvP%s", d1, cc[k]), &p1[k]);
        t->SetBranchStatus(Form("rc%sSvP%s", d2, cc[k]), 1); t->SetBranchAddress(Form("rc%sSvP%s", d2, cc[k]), &p2[k]);
      }
      for (auto& b : std::vector<std::pair<std::string, float*>>{{"wLevy", &w}, {"rcCosTheta", &cosT}, {Form("rc%sDcaXY", d1), &dx1}, {Form("rc%sDcaXY", d2), &dx2}}) {
        t->SetBranchStatus(b.first.c_str(), 1); t->SetBranchAddress(b.first.c_str(), b.second);
      }
      t->SetBranchStatus("passAll", 1); t->SetBranchAddress("passAll", &pass);
      TRandom3 r(7);
      const double m1 = mPi, m2 = ks ? mPi : mP;
      const Long64_t n = std::min<Long64_t>(t->GetEntries(), 3000000);
      for (Long64_t i = 0; i < n; ++i) {
        t->GetEntry(i);
        if (!pass) continue;
        double a1[3] = {p1[0], p1[1], p1[2]}, a2[3] = {p2[0], p2[1], p2[2]}, o1[3], o2[3], pt;
        smear(a1, m1, 0.0375, 0.014, 0.991, r.Gaus(), o1); smear(a2, m2, 0.0375, 0.014, 0.991, r.Gaus(), o2);
        const double m = minv(o1, m1, o2, m2, pt);
        const int k = bin(pt); if (k < 0) continue;
        const bool tight = cosT > diraMin && std::fabs(dx1) > piDca && std::fabs(dx2) > (ks ? piDca : pDca);
        for (int tt = 0; tt < (tight ? 2 : 1); ++tt) { hs[tt][k]->Fill(m, w); hs[tt][nb]->Fill(m, w); }
      }
    }
    for (int j = 0; j <= nb; ++j) {
      G db = gfit(hd[0][j], ks), dt = gfit(hd[1][j], ks), sb = gfit(hs[0][j], ks), st = gfit(hs[1][j], ks);
      const double rd = dt.s / db.s, rs = st.s / sb.s;
            const double fr = hd[1][j]->GetEntries() / std::max(1., hd[0][j]->GetEntries());
      const double erd = rd * std::sqrt(std::pow(dt.es / dt.s, 2) + std::pow(db.es / db.s, 2));   // conservative (ignores overlap)
      printf("BIAS %s pt %.1f-%.1f data: sigma base %.2f+-%.2f tight %.2f+-%.2f ratio %.3f+-%.3f (kept %.2f of entries) | DDSim: base %.2f tight %.2f ratio %.3f | mean data base %.2f tight %.2f\n",
             ks ? "Ks" : "Lambda", j < nb ? pb[j] : pb[0], j < nb ? pb[j + 1] : pb[nb], db.s, db.es, dt.s, dt.es, rd, erd, fr, sb.s, st.s, rs, db.m, dt.m);
    }
  }
}
