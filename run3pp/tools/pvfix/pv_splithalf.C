// pv_splithalf.C : PV resolution of the pair-average method (PHSimpleVertexFinder-like) and the KFParticle vertex fit,
// from split-half vertices: the tracks of each SvtxVertex (same quality cuts as the d0-pvfix refits) are split by
// index parity, each half is fitted with each method, and d = v(half1) - v(half2) is histogrammed. For a method,
// sigma(d)/2 estimates the resolution of the full vertex; for the KF fit d/sqrt(C1+C2) tests its error scale.
//   root -l -b -q 'pv_splithalf.C+("<dst list>", <max events per file>, <min tracks>)'
#define HomogeneousField  // as in the KFParticle_sPHENIX build (-DHomogeneousField)
#include <globalvertex/SvtxVertex.h>
#include <globalvertex/SvtxVertexMap.h>
#include <trackbase/TrkrDefs.h>
#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/TrackSeed.h>
#include <trackbase_historic/TrackSeedContainer.h>

#include <KFParticle.h>
#include <KFVertex.h>

#include <TFile.h>
#include <TTree.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <memory>
#include <set>
#include <string>
#include <vector>

R__LOAD_LIBRARY(libtrackbase_historic_io.so)
R__LOAD_LIBRARY(libglobalvertex_io.so)
R__LOAD_LIBRARY(libKFParticle.so)

namespace
{
  struct Trk
  {
    SvtxTrack *t;
    double a[3], b[3], cov[3][3];
  };
  double dot(const double *u, const double *v) { return u[0] * v[0] + u[1] * v[1] + u[2] * v[2]; }
  double dcaLines(const double *a1, const double *b1, const double *a2, const double *b2, double *p1, double *p2)
  {
    double bx[3] = {b1[1] * b2[2] - b1[2] * b2[1], b1[2] * b2[0] - b1[0] * b2[2], b1[0] * b2[1] - b1[1] * b2[0]};
    double mag = std::sqrt(dot(bx, bx));
    if (mag == 0) return 999;
    double am[3] = {a2[0] - a1[0], a2[1] - a1[1], a2[2] - a1[2]};
    double dca = dot(bx, am) / mag, b12 = dot(b2, b1);
    double X = dot(b1, b2) - dot(b1, b1) * dot(b2, b2) / b12;
    double Y = (dot(a2, b2) - dot(a1, b2)) - (dot(a2, b1) - dot(a1, b1)) * dot(b2, b2) / b12;
    double c = Y / X, d = c * dot(b1, b1) / b12 - (dot(a2, b1) - dot(a1, b1)) / b12;
    for (int k = 0; k < 3; ++k)
    {
      p1[k] = a1[k] + c * b1[k];
      p2[k] = a2[k] + d * b2[k];
    }
    return dca;
  }
  double median(std::vector<double> v)
  {
    std::sort(v.begin(), v.end());
    size_t n = v.size();
    return n % 2 ? v[n / 2] : 0.5 * (v[n / 2] + v[n / 2 - 1]);
  }
  bool pairVertex(const std::vector<Trk> &t, double *pos)
  {
    std::vector<std::array<double, 6>> pc;
    for (double cut : {0.05, 0.15})
    {
      pc.clear();
      for (size_t i = 0; i < t.size(); ++i)
        for (size_t j = i + 1; j < t.size(); ++j)
        {
          double p1[3], p2[3];
          double d = dcaLines(t[i].a, t[i].b, t[j].a, t[j].b, p1, p2);
          if (std::fabs(d) < cut && std::fabs(p1[0]) < 0.2 && std::fabs(p1[1]) < 0.2 && std::fabs(p2[0]) < 0.2 && std::fabs(p2[1]) < 0.2)
            pc.push_back({p1[0], p1[1], p1[2], p2[0], p2[1], p2[2]});
        }
      if (!pc.empty()) break;
    }
    if (pc.empty()) return false;
    std::vector<double> vx, vy, vz;
    for (auto &p : pc)
    {
      vx.insert(vx.end(), {p[0], p[3]});
      vy.insert(vy.end(), {p[1], p[4]});
      vz.insert(vz.end(), {p[2], p[5]});
    }
    double mx = median(vx), my = median(vy), mz = median(vz), s[3] = {0, 0, 0}, w = 0;
    if (vx.size() < 3)
    {
      pos[0] = 0.5 * (vx[0] + vx[1]);
      pos[1] = 0.5 * (vy[0] + vy[1]);
      pos[2] = 0.5 * (vz[0] + vz[1]);
      return true;
    }
    for (auto &p : pc)
      if (std::fabs(p[0] - mx) < 0.1 && std::fabs(p[1] - my) < 0.1 && std::fabs(p[3] - mx) < 0.1 && std::fabs(p[4] - my) < 0.1)
      {
        for (int k = 0; k < 3; ++k) s[k] += p[k] + p[k + 3];
        w += 2;
      }
    if (w > 0)
      for (int k = 0; k < 3; ++k) pos[k] = s[k] / w;
    else
    {
      pos[0] = mx;
      pos[1] = my;
      pos[2] = mz;
    }
    return true;
  }
  float gErrScale = 1;  // track covariance scaled by gErrScale^2
  bool kfVertex(const std::vector<Trk> &t, float chiCut, double *pos, double *var, double &chi2ndf, int &nUsed)
  {
    size_t n = t.size();
    std::vector<KFParticle> parts(n);
    std::vector<const KFParticle *> ptr(n);
    std::unique_ptr<Bool_t[]> flag(new Bool_t[n]);
    for (size_t i = 0; i < n; ++i)
    {
      const SvtxTrack *s = t[i].t;
      float par[6] = {s->get_x(), s->get_y(), s->get_z(), s->get_px(), s->get_py(), s->get_pz()}, cov[21];
      int k = 0;
      for (int a = 0; a < 6; ++a)
        for (int b = 0; b <= a; ++b) cov[k++] = s->get_error(a, b) * gErrScale * gErrScale;
      parts[i].Create(par, cov, (Int_t) s->get_charge(), -1);
      ptr[i] = &parts[i];
    }
    KFVertex v;
    v.ConstructPrimaryVertex(ptr.data(), (int) n, flag.get(), chiCut);
    nUsed = 0;
    for (size_t i = 0; i < n; ++i) nUsed += flag[i] ? 1 : 0;
    if (nUsed < 2 || !std::isfinite(v.GetX())) return false;
    pos[0] = v.GetX();
    pos[1] = v.GetY();
    pos[2] = v.GetZ();
    var[0] = v.GetCovariance(0, 0);
    var[1] = v.GetCovariance(1, 1);
    var[2] = v.GetCovariance(2, 2);
    chi2ndf = v.GetNDF() > 0 ? v.GetChi2() / v.GetNDF() : -1;
    return true;
  }
  double robust(std::vector<double> v)  // 1.4826 * MAD
  {
    if (v.size() < 5) return NAN;
    double m = median(v);
    for (auto &x : v) x = std::fabs(x - m);
    return 1.4826 * median(v);
  }
}  // namespace

void pv_splithalf(const char *list, int maxEv = 200, int minTrk = 6, float errScale = 1)
{
  gErrScale = errScale;
  printf("track error scale %.2f\n", errScale);
  KFParticle::SetField(14.0);  // kG, sPHENIX solenoid at the origin (KFParticle_sPHENIX reads it from the field map)
  std::ifstream in(list);
  std::string fn;
  const std::vector<std::string> meth = {"pair", "KF chi2cut 1e9", "KF chi2cut 9", "KF chi2cut 3.5"};
  const float cuts[4] = {0, 1e9, 9, 3.5};
  std::vector<std::vector<double>> dx(4), dz(4), px(4), pz(4), fx(4), fz(4), chi(4), used(4);
  long nv = 0, nvSel = 0;
  while (std::getline(in, fn))
  {
    TFile *f = TFile::Open(fn.c_str());
    if (!f || f->IsZombie()) continue;
    TTree *T = (TTree *) f->Get("T");
    SvtxTrackMap *tm = nullptr;
    SvtxVertexMap *vm = nullptr;
    TrackSeedContainer *sic = nullptr;
    T->SetBranchAddress("DST#SVTX#SvtxTrackMap", &tm);
    T->SetBranchAddress("DST#SVTX#SvtxVertexMap", &vm);
    T->SetBranchAddress("DST#SVTX#SiliconTrackSeedContainer", &sic);
    for (long ie = 0; ie < std::min<long>(T->GetEntries(), maxEv); ++ie)
    {
      T->GetEntry(ie);
      for (auto &[vk, v] : *vm)
      {
        ++nv;
        std::vector<Trk> trks;
        for (auto it = v->begin_tracks(); it != v->end_tracks(); ++it)
        {
          SvtxTrack *t = tm->get(*it);
          if (!t || t->get_quality() > 300 || t->get_pt() < 0.1) continue;
          const TrackSeed *ss = (sic && t->get_silicon_seed_id() != USHRT_MAX) ? sic->get(t->get_silicon_seed_id()) : t->get_silicon_seed();
          int nm = 0;
          if (ss)
            for (auto k = ss->begin_cluster_keys(); k != ss->end_cluster_keys(); ++k) nm += TrkrDefs::getTrkrId(*k) == TrkrDefs::mvtxId;
          if (nm < 3) continue;
          Trk r{t, {t->get_x(), t->get_y(), t->get_z()}, {}, {}};
          double p = t->get_p();
          r.b[0] = t->get_px() / p;
          r.b[1] = t->get_py() / p;
          r.b[2] = t->get_pz() / p;
          trks.push_back(r);
        }
        if ((int) trks.size() < minTrk) continue;
        ++nvSel;
        std::vector<Trk> h[2];
        for (size_t i = 0; i < trks.size(); ++i) h[i % 2].push_back(trks[i]);
        for (int m = 0; m < 4; ++m)
        {
          double p0[3], p1[3], pf[3], v0[3] = {0, 0, 0}, v1[3] = {0, 0, 0}, vf[3], c2 = -1, c0, c1;
          int u0 = 0, u1 = 0, uf = 0;
          bool ok0, ok1, okf;
          if (m == 0)
          {
            ok0 = pairVertex(h[0], p0);
            ok1 = pairVertex(h[1], p1);
            okf = pairVertex(trks, pf);
          }
          else
          {
            ok0 = kfVertex(h[0], cuts[m], p0, v0, c0, u0);
            ok1 = kfVertex(h[1], cuts[m], p1, v1, c1, u1);
            okf = kfVertex(trks, cuts[m], pf, vf, c2, uf);
          }
          if (okf)
          {
            fx[m].push_back(1e4 * (pf[0] - v->get_x()));
            fz[m].push_back(1e4 * (pf[2] - v->get_z()));
            if (m) chi[m].push_back(c2), used[m].push_back((double) uf / trks.size());
          }
          if (!ok0 || !ok1) continue;
          dx[m].push_back(1e4 * (p0[0] - p1[0]));
          dz[m].push_back(1e4 * (p0[2] - p1[2]));
          if (m)
          {
            px[m].push_back((p0[0] - p1[0]) / std::sqrt(v0[0] + v1[0]));
            pz[m].push_back((p0[2] - p1[2]) / std::sqrt(v0[2] + v1[2]));
          }
        }
      }
    }
    f->Close();
  }
  printf("vertices %ld, with >= %d selected tracks %ld\n", nv, minTrk, nvSel);
  printf("%-16s %8s | split-half: sigma(dx)/2 sigma(dz)/2 [um, robust] | KF pull width x z | full - stored: sigma x z [um] | KF chi2/ndf median, track fraction used\n", "method", "N");
  for (int m = 0; m < 4; ++m)
    printf("%-16s %8zu | %8.1f %8.1f | %6.2f %6.2f | %8.1f %8.1f | %6.2f %5.2f\n", meth[m].c_str(), dx[m].size(), robust(dx[m]) / 2, robust(dz[m]) / 2,
           m ? robust(px[m]) : NAN, m ? robust(pz[m]) : NAN, robust(fx[m]), robust(fz[m]), m ? median(chi[m]) : NAN, m ? median(used[m]) : NAN);
}
