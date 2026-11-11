// kc_stats.C : summary of the kappa correction on the kc skim (kappa_apply.C output): kappa per charge, mass and pT
// changes (quantile-based: median and the half-width of the central 68% interval, plus the fraction of candidates beyond
// 20 MeV, because a few far outliers dominate a plain r.m.s.), candidate-set check. Rows: KC <key> <value>
//   root -l -b -q 'kc_stats.C+("<kc skim>")'
#include <ROOT/RDataFrame.hxx>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <vector>

namespace
{
  void qstats(const char *key, std::vector<float> v, double tail)
  {
    std::sort(v.begin(), v.end());
    auto at = [&](double p) { return v[std::min(v.size() - 1, (size_t) (p * v.size()))]; };
    long nt = std::count_if(v.begin(), v.end(), [&](float x) { return std::fabs(x) > tail; });
    printf("KC %sMed %.2f\nKC %sHw %.2f\nKC %sTail %.4f\n", key, at(0.5), key, 0.5 * (at(0.84) - at(0.16)), key, (double) nt / v.size());
  }
}  // namespace

void kc_stats(const char *in)
{
  ROOT::EnableImplicitMT(6);
  ROOT::RDataFrame d("DecayTree", in);
  auto n = d.Count();
  auto nfin = d.Filter("std::isfinite(D0_mass_kc) && std::isfinite(D0_mass_rc)").Count();
  auto kp = d.Filter("track_1_charge>0").Stats("track_1_kappa"), km = d.Filter("track_1_charge<0").Stats("track_1_kappa");
  auto rc = d.Filter("D0_decayLength<0.1").Define("x", "1e3f*(D0_mass_rc-D0_mass)").Take<float>("x");
  auto kc = d.Define("x", "1e3f*(D0_mass_kc-D0_mass_rc)").Take<float>("x");
  auto pt = d.Filter("D0_pT>0").Define("x", "D0_pT_kc/D0_pT-1.f").Take<float>("x");
  auto npt0 = d.Filter("!(D0_pT>0)").Count();
  printf("KC N %llu\nKC Nfinite %llu\nKC NptZero %llu\n", *n, *nfin, *npt0);
  printf("KC kappaPlusMean %.4f\nKC kappaPlusRms %.4f\nKC kappaMinusMean %.4f\nKC kappaMinusRms %.4f\n", kp->GetMean(), kp->GetRMS(), km->GetMean(), km->GetRMS());
  qstats("dmRc", *rc, 20);   // recomputed - stored mass, L < 1 mm, MeV
  qstats("dmKc", *kc, 20);   // corrected - recomputed mass, MeV
  qstats("ptRel", *pt, 0.1); // relative pair-pT change
}
