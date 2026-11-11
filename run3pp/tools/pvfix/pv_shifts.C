// pv_shifts.C : primary-vertex versions in the d0-pvfix v2 ntuples (one row per candidate).
//  (a) closure of the all-track refits vs the stored PV, in bins of the number of tracks of the refit:
//      pair refit (refitAll) and KF fit (kfAll): fraction with |dxy| < 1 um, dxy quantiles, robust dz width
//  (b) effect of removing the daughters (fix 2/3): refitAll -> refit and kfAll -> kf shift, split by how many daughters
//      are attached to the stored PV
//  (c) KF fit status/track usage
// Prints machine-readable rows: CLOS <method> <nbin> N f1um q50 q90 dzw ; SHIFT <method> <attached> N q50 q90 dzw
//   root -l -b -q 'pv_shifts.C+("<file glob>")'
#include <ROOT/RDataFrame.hxx>
#include <TH1D.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

namespace
{
  double q(std::vector<float> v, double p)
  {
    if (v.empty()) return NAN;
    std::sort(v.begin(), v.end());
    return v[std::min(v.size() - 1, (size_t) (p * v.size()))];
  }
  double robust(std::vector<float> v)
  {
    if (v.size() < 5) return NAN;
    double m = q(v, 0.5);
    for (auto &x : v) x = std::fabs(x - m);
    return 1.4826 * q(v, 0.5);
  }
}  // namespace

void pv_shifts(const char *in)
{
  ROOT::EnableImplicitMT(6);
  ROOT::RDataFrame d("DecayTree", in);
  printf("candidates %llu\n", *d.Count());
  struct M
  {
    std::string name, pre, nt, st;
  };
  for (M m : {M{"pair", "primary_vertex_refitAll", "primary_vertex_refitAll_nTracks", "primary_vertex_refitAll_status"},
              M{"kf", "primary_vertex_kfAll", "primary_vertex_kfAll_nTracks", "primary_vertex_kfAll_status"}})
  {
    auto x = d.Filter(m.st + "==1")
                 .Define("dxy", "1e4f*(float)std::hypot(" + m.pre + "_x-primary_vertex_x," + m.pre + "_y-primary_vertex_y)")
                 .Define("dz", "1e4f*(" + m.pre + "_z-primary_vertex_z)")
                 .Define("nt", "(int)" + m.nt);
    const char *bl[4] = {"2", "3", "4-6", ">=7"};
    const char *bc[4] = {"nt==2", "nt==3", "nt>=4 && nt<=6", "nt>=7"};
    for (int b = 0; b < 4; ++b)
    {
      auto y = x.Filter(bc[b]);
      auto vxy = y.Take<float>("dxy"), vz = y.Take<float>("dz");
      std::vector<float> a = *vxy, z = *vz;
      long n1 = std::count_if(a.begin(), a.end(), [](float v) { return v < 1; });
      printf("CLOS %s %s N %zu f1um %.3f q50 %.1f q90 %.1f dzw %.1f\n", m.name.c_str(), bl[b], a.size(), a.empty() ? 0. : (double) n1 / a.size(), q(a, 0.5), q(a, 0.9),
             robust(z));
    }
  }
  struct S
  {
    std::string name, all, no, st1, st2;
  };
  for (S s : {S{"pair", "primary_vertex_refitAll", "primary_vertex_refit", "primary_vertex_refitAll_status", "primary_vertex_refit_status"},
              S{"kf", "primary_vertex_kfAll", "primary_vertex_kf", "primary_vertex_kfAll_status", "primary_vertex_kf_status"}})
  {
    auto x = d.Filter(s.st1 + "==1 && " + s.st2 + "==1")
                 .Define("sxy", "1e4f*(float)std::hypot(" + s.no + "_x-" + s.all + "_x," + s.no + "_y-" + s.all + "_y)")
                 .Define("sz", "1e4f*(" + s.no + "_z-" + s.all + "_z)")
                 .Define("att", "(int)(track_1_inPV==1)+(int)(track_2_inPV==1)");
    for (int a = 2; a >= 0; --a)
    {
      auto y = x.Filter(Form("att==%d", a));
      std::vector<float> v = *y.Take<float>("sxy"), z = *y.Take<float>("sz");
      printf("SHIFT %s %d N %zu q50 %.1f q90 %.1f dzw %.1f\n", s.name.c_str(), a, v.size(), q(v, 0.5), q(v, 0.9), robust(z));
    }
  }
  auto nAll = d.Count();
  auto okA = d.Filter("primary_vertex_kfAll_status==1").Count(), ok = d.Filter("primary_vertex_kf_status==1").Count();
  auto okR = d.Filter("primary_vertex_refit_status==1").Count(), okRA = d.Filter("primary_vertex_refitAll_status==1").Count();
  auto use = d.Filter("primary_vertex_kfAll_status==1 && primary_vertex_refitAll_status==1 && primary_vertex_refitAll_nTracks>0")
                 .Define("u", "(float)primary_vertex_kfAll_nTracks/(float)primary_vertex_refitAll_nTracks")
                 .Mean("u");
  printf("STATUS N %llu refitAll %.4f refit %.4f kfAll %.4f kf %.4f kfUse %.3f\n", *nAll, (double) *okRA / *nAll, (double) *okR / *nAll, (double) *okA / *nAll,
         (double) *ok / *nAll, *use);
}
