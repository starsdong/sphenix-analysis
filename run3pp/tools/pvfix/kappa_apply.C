// kappa_apply.C : apply the momentum-scale map (kappa_lookup.csv, MomentumScaleMap/run79509to79516_20261102) to a D0
// skim and write a copy with additional branches; all original branches are kept unchanged.
//   pT_corr = kappa(q, pT, eta, phi) * pT  (direction unchanged, so p, px, py, pz scale by kappa as well)
//   kappa: trilinear interpolation in (x = ln pT, eta, phi) on the map grid, separately for q = -1 / +1;
//          pT clamped to [0.1, 20] GeV, eta to [-1.5, 1.5]; phi periodic.
// New branches:
//   track_{1,2}_kappa, track_{1,2}_pT_kc
//   D0_mass_rc : invariant mass recomputed from the stored daughter (pT, eta, phi) with the stored mass hypotheses,
//                no correction (reference: isolates the effect of the map from the recomputation)
//   D0_mass_kc, D0_pT_kc : the same with corrected daughter momenta
//   D0_mass_kc_swap : corrected mass with the K/pi hypotheses exchanged
//   root -l -b -q 'kappa_apply.C+("<map csv>", "<in skim>", "<out file>")'
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>
#include <Math/Vector4D.h>
#include <TMath.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>

namespace
{
  struct KMap
  {
    std::vector<double> x, eta, phi;  // grid axes (x = ln pT)
    std::vector<float> k[2];          // [q index][ix][ieta][iphi]
    size_t idx(size_t i, size_t j, size_t l) const { return (i * eta.size() + j) * phi.size() + l; }
    static size_t lo(const std::vector<double> &a, double v, double &f)
    {
      if (v <= a.front()) { f = 0; return 0; }
      if (v >= a.back()) { f = 1; return a.size() - 2; }
      size_t i = std::upper_bound(a.begin(), a.end(), v) - a.begin() - 1;
      f = (v - a[i]) / (a[i + 1] - a[i]);
      return i;
    }
    float eval(int q, float pt, float et, float ph) const
    {
      const auto &K = k[q > 0 ? 1 : 0];
      double fx, fe, fp;
      double xv = std::log(std::min(std::max((double) pt, 0.1), 20.0));
      double ev = std::min(std::max((double) et, -1.5), 1.5);
      double pv = std::remainder((double) ph, 2 * M_PI);  // [-pi, pi]
      size_t i = lo(x, xv, fx), j = lo(eta, ev, fe), l = lo(phi, pv, fp);
      double r = 0;
      for (int a = 0; a < 2; ++a)
        for (int b = 0; b < 2; ++b)
          for (int c = 0; c < 2; ++c)
            r += (a ? fx : 1 - fx) * (b ? fe : 1 - fe) * (c ? fp : 1 - fp) * K[idx(i + a, j + b, l + c)];
      return (float) r;
    }
  };
  KMap gMap;
  void loadMap(const char *csv)
  {
    std::ifstream in(csv);
    std::string line;
    std::vector<std::vector<double>> rows;
    std::set<double> sx, se, sp;
    while (std::getline(in, line))
    {
      if (line.empty() || line[0] == '#' || line[0] == 'q') continue;
      std::stringstream ss(line);
      std::vector<double> v;
      std::string tok;
      while (std::getline(ss, tok, ',')) v.push_back(std::stod(tok));
      // q,pT,curv,eta,phi,kappa,...
      rows.push_back(v);
      sx.insert(std::log(v[1]));
      se.insert(v[3]);
      sp.insert(v[4]);
    }
    auto round = [](std::set<double> &s) {  // merge values equal to 1e-6
      std::vector<double> o;
      for (double v : s)
        if (o.empty() || std::fabs(v - o.back()) > 1e-6) o.push_back(v);
      return o;
    };
    gMap.x = round(sx);
    gMap.eta = round(se);
    gMap.phi = round(sp);
    size_t n = gMap.x.size() * gMap.eta.size() * gMap.phi.size();
    for (auto &k : gMap.k) k.assign(n, NAN);
    auto find = [](const std::vector<double> &a, double v) {
      size_t i = std::lower_bound(a.begin(), a.end(), v - 1e-6) - a.begin();
      return i;
    };
    for (auto &v : rows)
    {
      size_t i = find(gMap.x, std::log(v[1])), j = find(gMap.eta, v[3]), l = find(gMap.phi, v[4]);
      gMap.k[v[0] > 0 ? 1 : 0][gMap.idx(i, j, l)] = (float) v[5];
    }
    long nan = 0;
    for (auto &k : gMap.k)
      for (float f : k) nan += std::isnan(f);
    printf("kappa map: %zu rows, grid ln(pT) %zu [%.3f, %.3f], eta %zu [%.2f, %.2f], phi %zu [%.4f, %.4f], unfilled nodes %ld\n", rows.size(), gMap.x.size(),
           gMap.x.front(), gMap.x.back(), gMap.eta.size(), gMap.eta.front(), gMap.eta.back(), gMap.phi.size(), gMap.phi.front(), gMap.phi.back(), nan);
  }
  float kappaOf(int q, float pt, float eta, float phi) { return gMap.eval(q, pt, eta, phi); }
  float massOf(float pt1, float eta1, float phi1, float m1, float pt2, float eta2, float phi2, float m2)
  {
    ROOT::Math::PtEtaPhiMVector a(pt1, eta1, phi1, m1), b(pt2, eta2, phi2, m2);
    return (float) (a + b).M();
  }
  float pairPt(float pt1, float phi1, float pt2, float phi2)
  {
    return (float) std::hypot(pt1 * std::cos(phi1) + pt2 * std::cos(phi2), pt1 * std::sin(phi1) + pt2 * std::sin(phi2));
  }
  const float mK = 0.493677, mPi = 0.13957;
  float hypo(int pdg) { return std::abs(pdg) == 321 ? mK : mPi; }
}  // namespace

void kappa_apply(const char *csv, const char *in, const char *out)
{
  loadMap(csv);
  // spot checks against grid nodes
  printf("spot check kappa(q=-1, pT=2.1, eta=0, phi=0) = %.5f, kappa(q=+1, 2.1, 0, 0) = %.5f\n", kappaOf(-1, 2.1, 0, 0), kappaOf(1, 2.1, 0, 0));
  ROOT::EnableImplicitMT(6);
  ROOT::RDataFrame d("DecayTree", in);
  auto o = d.Define("track_1_kappa", [](char c, float p, float e, float f) { return kappaOf(c > 0 ? 1 : -1, p, e, f); },
                    {"track_1_charge", "track_1_pT", "track_1_pseudorapidity", "track_1_phi"})
               .Define("track_2_kappa", [](char c, float p, float e, float f) { return kappaOf(c > 0 ? 1 : -1, p, e, f); },
                       {"track_2_charge", "track_2_pT", "track_2_pseudorapidity", "track_2_phi"})
               .Define("track_1_pT_kc", "track_1_kappa*track_1_pT")
               .Define("track_2_pT_kc", "track_2_kappa*track_2_pT")
               .Define("m1h", [](int p) { return hypo(p); }, {"track_1_PDG_ID"})
               .Define("m2h", [](int p) { return hypo(p); }, {"track_2_PDG_ID"})
               .Define("D0_mass_rc", massOf, {"track_1_pT", "track_1_pseudorapidity", "track_1_phi", "m1h", "track_2_pT", "track_2_pseudorapidity", "track_2_phi", "m2h"})
               .Define("D0_mass_kc", massOf,
                       {"track_1_pT_kc", "track_1_pseudorapidity", "track_1_phi", "m1h", "track_2_pT_kc", "track_2_pseudorapidity", "track_2_phi", "m2h"})
               .Define("D0_mass_kc_swap", massOf,
                       {"track_1_pT_kc", "track_1_pseudorapidity", "track_1_phi", "m2h", "track_2_pT_kc", "track_2_pseudorapidity", "track_2_phi", "m1h"})
               .Define("D0_pT_kc", pairPt, {"track_1_pT_kc", "track_1_phi", "track_2_pT_kc", "track_2_phi"});
  std::vector<std::string> cols;
  for (auto &c : d.GetColumnNames()) cols.push_back(c);
  for (const char *c : {"track_1_kappa", "track_2_kappa", "track_1_pT_kc", "track_2_pT_kc", "D0_mass_rc", "D0_mass_kc", "D0_mass_kc_swap", "D0_pT_kc"}) cols.push_back(c);
  auto hk = o.Histo1D({"hk", "", 200, 0.6, 1.4}, "track_1_kappa");
  auto hdr = o.Define("dm", "1e3*(D0_mass_rc-D0_mass)").Filter("D0_decayLength<0.1").Histo1D({"hdm", "", 400, -100, 100}, "dm");
  auto hdk = o.Define("dk", "1e3*(D0_mass_kc-D0_mass_rc)").Histo1D({"hdk", "", 400, -200, 200}, "dk");
  ROOT::RDF::RSnapshotOptions so;
  so.fCompressionAlgorithm = ROOT::kZSTD;
  so.fCompressionLevel = 5;
  so.fLazy = true;
  auto snap = o.Snapshot("DecayTree", out, cols, so);
  *snap;
  printf("kappa (track 1): mean %.4f rms %.4f\n", hk->GetMean(), hk->GetRMS());
  printf("D0_mass_rc - D0_mass (L < 1 mm): mean %.2f MeV, rms %.2f MeV\n", hdr->GetMean(), hdr->GetRMS());
  printf("D0_mass_kc - D0_mass_rc: mean %.2f MeV, rms %.2f MeV\n", hdk->GetMean(), hdk->GetRMS());
  printf("kappa_apply: wrote %s with %zu branches\n", out, cols.size());
}
