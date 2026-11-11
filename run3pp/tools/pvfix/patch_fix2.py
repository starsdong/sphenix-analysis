#!/usr/bin/env python3
"""d0-pvfix fix 2: refit the candidate's primary vertex without its daughter tracks and store the PV-dependent
quantities w.r.t. the refitted PV in new *_PVrefit branches (original branches unchanged).
The refit re-runs the PHSimpleVertexFinder algorithm (ana.581 production settings, rootmacros/Trkr_Reco.C:
pair DCA < 0.05 cm (x3 if no pair), pT > 0.1, quality < 300, >= 3 MVTX clusters, pair PCA |x|,|y| < 0.2 cm,
outlier pair cut 0.10 cm around the median, position = mean of the accepted pair PCA points, covariance =
sum of track position covariances / sqrt(n)) on the tracks attached to the SvtxVertex:
  primary_vertex_refitAll_*  all attached tracks (closure test against the stored vertex)
  primary_vertex_refit_*     attached tracks minus the candidate's daughters (used for the *_PVrefit quantities)
  *_PVrefitAll quantities    w.r.t. the all-tracks refit: original -> PVrefitAll = algorithm difference,
                             PVrefitAll -> PVrefit = effect of excluding the daughters
KFP_PVREFIT_CLOSURE=1 in the environment makes the *_PVrefit quantities use the stored PV (code-path closure).
Run in offline/packages/KFParticle_sPHENIX."""

def rd(p):
    b = open(p, 'rb').read().decode()
    return b, ('\r\n' if '\r\n' in b else '\n')

def wr(p, s):
    open(p, 'wb').write(s.encode())

def sub(s, a, b, nl):
    a, b = a.replace('\n', nl), b.replace('\n', nl)
    assert s.count(a) == 1, ('not unique / missing', a[:80], s.count(a))
    return s.replace(a, b)

h, nl = rd('KFParticle_nTuple.h')
h = sub(h, '  int m_nTracksOfVertex{-1};\n', '''  int m_nTracksOfVertex{-1};

  // d0-pvfix fix 2: PV refit without the candidate's daughter tracks
  void pvRefit(PHCompositeNode* topNode, const KFParticle& mother, const KFParticle& pvIn, const std::vector<KFParticle>& daughters);
  static const int kMaxRefitTracks = 20;
  int m_nTracksOfVertex_noDaughters{-1};
  int m_pvr_status{-1};
  int m_pvr_ntrk{-1};
  int m_pvr_npair{-1};
  float m_pvr_x{0};
  float m_pvr_y{0};
  float m_pvr_z{0};
  float m_pvr_cov[6]{0};
  int m_pvra_status{-1};
  int m_pvra_ntrk{-1};
  int m_pvra_npair{-1};
  float m_pvra_x{0};
  float m_pvra_y{0};
  float m_pvra_z{0};
  float m_mr_dira{0};
  float m_mr_dira_xy{0};
  float m_mr_fdchi2{0};
  float m_mr_pvdca{0};
  float m_mr_pvdca_xy{0};
  float m_mr_L{0};
  float m_mr_Lerr{0};
  float m_mr_Lxy{0};
  float m_mr_Lxyerr{0};
  int m_dr_inPV[kMaxRefitTracks]{0};
  float m_dr_dca[kMaxRefitTracks]{0};
  float m_dr_dca_xy[kMaxRefitTracks]{0};
  float m_dr_dca_sig[kMaxRefitTracks]{0};
  float m_dr_dca_sig_xy[kMaxRefitTracks]{0};
  // same quantities w.r.t. the all-tracks refit (separates the algorithm difference from the daughter exclusion)
  float m_ma_dira{0};
  float m_ma_dira_xy{0};
  float m_ma_fdchi2{0};
  float m_ma_pvdca{0};
  float m_ma_pvdca_xy{0};
  float m_ma_L{0};
  float m_ma_Lerr{0};
  float m_ma_Lxy{0};
  float m_ma_Lxyerr{0};
  float m_da_dca[kMaxRefitTracks]{0};
  float m_da_dca_xy[kMaxRefitTracks]{0};
  float m_da_dca_sig[kMaxRefitTracks]{0};
  float m_da_dca_sig_xy[kMaxRefitTracks]{0};
''', nl)
wr('KFParticle_nTuple.h', h)

c, nl = rd('KFParticle_nTuple.cc')
c = sub(c, '#include <trackbase_historic/SvtxTrackMap.h>\n', '''#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/TrackSeed.h>           // d0-pvfix fix 2
#include <trackbase_historic/TrackSeedContainer.h>  // d0-pvfix fix 2
#include <trackbase/TrkrDefs.h>                     // d0-pvfix fix 2
#include <globalvertex/SvtxVertex.h>                // d0-pvfix fix 2
#include <globalvertex/SvtxVertexMap.h>             // d0-pvfix fix 2
#include <array>
#include <limits>
#include <set>
''', nl)

# branches
c = sub(c, '''    m_tree->Branch("primary_vertex_Covariance", &m_calculated_vertex_cov, "primary_vertex_Covariance[6]/F", 6);
''', '''    m_tree->Branch("primary_vertex_Covariance", &m_calculated_vertex_cov, "primary_vertex_Covariance[6]/F", 6);

    // d0-pvfix fix 2: PV refit without the candidate's daughters (see KFParticle_nTuple::pvRefit)
    auto bookF = [&](const std::string& n, float* a) { m_tree->Branch(n.c_str(), a, (n + "/F").c_str()); };
    auto bookI = [&](const std::string& n, int* a) { m_tree->Branch(n.c_str(), a, (n + "/I").c_str()); };
    bookF("primary_vertex_refit_x", &m_pvr_x);
    bookF("primary_vertex_refit_y", &m_pvr_y);
    bookF("primary_vertex_refit_z", &m_pvr_z);
    m_tree->Branch("primary_vertex_refit_Covariance", &m_pvr_cov, "primary_vertex_refit_Covariance[6]/F");
    bookI("primary_vertex_refit_status", &m_pvr_status);
    bookI("primary_vertex_refit_nTracks", &m_pvr_ntrk);
    bookI("primary_vertex_refit_nPairs", &m_pvr_npair);
    bookF("primary_vertex_refitAll_x", &m_pvra_x);
    bookF("primary_vertex_refitAll_y", &m_pvra_y);
    bookF("primary_vertex_refitAll_z", &m_pvra_z);
    bookI("primary_vertex_refitAll_status", &m_pvra_status);
    bookI("primary_vertex_refitAll_nTracks", &m_pvra_ntrk);
    bookI("primary_vertex_refitAll_nPairs", &m_pvra_npair);
    bookI("nTracksOfVertex_noDaughters", &m_nTracksOfVertex_noDaughters);
    bookF(mother_name + "_DIRA_PVrefit", &m_mr_dira);
    bookF(mother_name + "_DIRA_xy_PVrefit", &m_mr_dira_xy);
    bookF(mother_name + "_FDchi2_PVrefit", &m_mr_fdchi2);
    bookF(mother_name + "_PV_DCA_PVrefit", &m_mr_pvdca);
    bookF(mother_name + "_PV_DCA_xy_PVrefit", &m_mr_pvdca_xy);
    bookF(mother_name + "_decayLength_PVrefit", &m_mr_L);
    bookF(mother_name + "_decayLengthErr_PVrefit", &m_mr_Lerr);
    bookF(mother_name + "_decayLength_xy_PVrefit", &m_mr_Lxy);
    bookF(mother_name + "_decayLengthErr_xy_PVrefit", &m_mr_Lxyerr);
    bookF(mother_name + "_DIRA_PVrefitAll", &m_ma_dira);
    bookF(mother_name + "_DIRA_xy_PVrefitAll", &m_ma_dira_xy);
    bookF(mother_name + "_FDchi2_PVrefitAll", &m_ma_fdchi2);
    bookF(mother_name + "_PV_DCA_PVrefitAll", &m_ma_pvdca);
    bookF(mother_name + "_PV_DCA_xy_PVrefitAll", &m_ma_pvdca_xy);
    bookF(mother_name + "_decayLength_PVrefitAll", &m_ma_L);
    bookF(mother_name + "_decayLengthErr_PVrefitAll", &m_ma_Lerr);
    bookF(mother_name + "_decayLength_xy_PVrefitAll", &m_ma_Lxy);
    bookF(mother_name + "_decayLengthErr_xy_PVrefitAll", &m_ma_Lxyerr);
    for (int i = 0; i < std::min(m_num_tracks_nTuple, kMaxRefitTracks); ++i)
    {
      std::string dn = "track_" + std::to_string(i + 1);
      bookI(dn + "_inPV", &m_dr_inPV[i]);
      bookF(dn + "_PV_DCA_PVrefit", &m_dr_dca[i]);
      bookF(dn + "_PV_DCA_xy_PVrefit", &m_dr_dca_xy[i]);
      bookF(dn + "_PV_DCA_sig_PVrefit", &m_dr_dca_sig[i]);
      bookF(dn + "_PV_DCA_sig_xy_PVrefit", &m_dr_dca_sig_xy[i]);
      bookF(dn + "_PV_DCA_PVrefitAll", &m_da_dca[i]);
      bookF(dn + "_PV_DCA_xy_PVrefitAll", &m_da_dca_xy[i]);
      bookF(dn + "_PV_DCA_sig_PVrefitAll", &m_da_dca_sig[i]);
      bookF(dn + "_PV_DCA_sig_xy_PVrefitAll", &m_da_dca_sig_xy[i]);
    }
''', nl)

# pristine copy of the mother at the start of fillBranch
c = sub(c, '''  const float speedOfLight = 2.99792458e-2;

  KFParticle temp;
''', '''  const float speedOfLight = 2.99792458e-2;
  const KFParticle motherForRefit = motherParticle;  // d0-pvfix fix 2

  KFParticle temp;
''', nl)

c = sub(c, '''  else
  {
    m_nTracksOfVertex = 0;
  }
''', '''  else
  {
    m_nTracksOfVertex = 0;
  }
  if (m_constrain_to_vertex_nTuple && !m_use_fake_pv_nTuple)
  {
    pvRefit(topNode, motherForRefit, vertex_fillbranch, daughters);  // d0-pvfix fix 2
  }
''', nl)

c += '''
// ---------------------------------------------------------------- d0-pvfix fix 2
namespace
{
  struct PVfixTrack
  {
    unsigned int id;
    double a[3];
    double b[3];
    double cov[3][3];
  };

  double pvfixDot(const double* u, const double* v) { return u[0] * v[0] + u[1] * v[1] + u[2] * v[2]; }

  // copy of PHSimpleVertexFinder::dcaTwoLines (straight lines a + c b)
  double pvfixDcaTwoLines(const double* a1, const double* b1, const double* a2, const double* b2, double* p1, double* p2)
  {
    double bx[3] = {b1[1] * b2[2] - b1[2] * b2[1], b1[2] * b2[0] - b1[0] * b2[2], b1[0] * b2[1] - b1[1] * b2[0]};
    double mag = std::sqrt(pvfixDot(bx, bx));
    if (mag == 0)
    {
      return 999;
    }
    double am[3] = {a2[0] - a1[0], a2[1] - a1[1], a2[2] - a1[2]};
    double dca = pvfixDot(bx, am) / mag;
    double b12 = pvfixDot(b2, b1);
    double X = pvfixDot(b1, b2) - pvfixDot(b1, b1) * pvfixDot(b2, b2) / b12;
    double Y = (pvfixDot(a2, b2) - pvfixDot(a1, b2)) - (pvfixDot(a2, b1) - pvfixDot(a1, b1)) * pvfixDot(b2, b2) / b12;
    double cc = Y / X;
    double F = pvfixDot(b1, b1) / b12;
    double G = -(pvfixDot(a2, b1) - pvfixDot(a1, b1)) / b12;
    double dd = cc * F + G;
    for (int k = 0; k < 3; ++k)
    {
      p1[k] = a1[k] + cc * b1[k];
      p2[k] = a2[k] + dd * b2[k];
    }
    return dca;
  }

  double pvfixMedian(std::vector<double> v)
  {
    std::sort(v.begin(), v.end());
    size_t n = v.size();
    return (n % 2 == 0) ? 0.5 * (v[n / 2] + v[n / 2 - 1]) : v[n / 2];
  }

  // PHSimpleVertexFinder position and covariance for one set of tracks; returns 1 if a position was found
  int pvfixSimpleVertex(const std::vector<PVfixTrack>& t, double* pos, double* cov, int& npair, int& ntrk)
  {
    const double beamlineCut = 0.2, outlierCut = 0.10;
    std::vector<std::array<double, 6>> pc;
    std::set<size_t> used;
    for (double cut : {0.05, 0.15})
    {
      pc.clear();
      used.clear();
      for (size_t i = 0; i < t.size(); ++i)
      {
        for (size_t j = i + 1; j < t.size(); ++j)
        {
          double p1[3], p2[3];
          double d = pvfixDcaTwoLines(t[i].a, t[i].b, t[j].a, t[j].b, p1, p2);
          if (std::fabs(d) < cut && std::fabs(p1[0]) < beamlineCut && std::fabs(p1[1]) < beamlineCut && std::fabs(p2[0]) < beamlineCut &&
              std::fabs(p2[1]) < beamlineCut)
          {
            pc.push_back({p1[0], p1[1], p1[2], p2[0], p2[1], p2[2]});
            used.insert(i);
            used.insert(j);
          }
        }
      }
      if (!pc.empty())
      {
        break;
      }
    }
    npair = pc.size();
    ntrk = used.size();
    if (pc.empty())
    {
      return 0;
    }
    std::vector<double> vx, vy, vz;
    for (auto& p : pc)
    {
      vx.insert(vx.end(), {p[0], p[3]});
      vy.insert(vy.end(), {p[1], p[4]});
      vz.insert(vz.end(), {p[2], p[5]});
    }
    if (vx.size() < 3)
    {
      pos[0] = 0.5 * (vx[0] + vx[1]);
      pos[1] = 0.5 * (vy[0] + vy[1]);
      pos[2] = 0.5 * (vz[0] + vz[1]);
    }
    else
    {
      double mx = pvfixMedian(vx), my = pvfixMedian(vy), mz = pvfixMedian(vz);
      double s[3] = {0, 0, 0};
      double w = 0;
      for (auto& p : pc)
      {
        if (std::fabs(p[0] - mx) < outlierCut && std::fabs(p[1] - my) < outlierCut && std::fabs(p[3] - mx) < outlierCut && std::fabs(p[4] - my) < outlierCut)
        {
          for (int k = 0; k < 3; ++k)
          {
            s[k] += p[k] + p[k + 3];
          }
          w += 2;
        }
      }
      if (w > 0)
      {
        for (int k = 0; k < 3; ++k)
        {
          pos[k] = s[k] / w;
        }
      }
      else
      {
        pos[0] = mx;
        pos[1] = my;
        pos[2] = mz;
      }
    }
    double c3[3][3] = {{0}};
    for (size_t i : used)
    {
      for (int a = 0; a < 3; ++a)
      {
        for (int b = 0; b < 3; ++b)
        {
          c3[a][b] += t[i].cov[a][b];
        }
      }
    }
    double sq = std::sqrt((double) used.size());
    int k = 0;
    for (int a = 0; a < 3; ++a)
    {
      for (int b = 0; b <= a; ++b)
      {
        cov[k++] = c3[a][b] / sq;
      }
    }
    return 1;
  }
}  // namespace

void KFParticle_nTuple::pvRefit(PHCompositeNode* topNode, const KFParticle& mother, const KFParticle& pvIn, const std::vector<KFParticle>& daughters)
{
  const float nan = std::numeric_limits<float>::quiet_NaN();
  m_nTracksOfVertex_noDaughters = -1;
  m_pvr_status = m_pvra_status = -1;
  m_pvr_ntrk = m_pvr_npair = m_pvra_ntrk = m_pvra_npair = -1;
  m_pvr_x = m_pvr_y = m_pvr_z = m_pvra_x = m_pvra_y = m_pvra_z = nan;
  for (float& x : m_pvr_cov)
  {
    x = nan;
  }
  m_mr_dira = m_mr_dira_xy = m_mr_fdchi2 = m_mr_pvdca = m_mr_pvdca_xy = m_mr_L = m_mr_Lerr = m_mr_Lxy = m_mr_Lxyerr = nan;
  m_ma_dira = m_ma_dira_xy = m_ma_fdchi2 = m_ma_pvdca = m_ma_pvdca_xy = m_ma_L = m_ma_Lerr = m_ma_Lxy = m_ma_Lxyerr = nan;
  const int nd = std::min((int) daughters.size(), kMaxRefitTracks);
  for (int i = 0; i < kMaxRefitTracks; ++i)
  {
    m_dr_inPV[i] = -1;
    m_dr_dca[i] = m_dr_dca_xy[i] = m_dr_dca_sig[i] = m_dr_dca_sig_xy[i] = nan;
    m_da_dca[i] = m_da_dca_xy[i] = m_da_dca_sig[i] = m_da_dca_sig_xy[i] = nan;
  }

  auto* vmap = findNode::getClass<SvtxVertexMap>(topNode, m_tools->get_vtx_map_node_name());
  auto* tmap = findNode::getClass<SvtxTrackMap>(topNode, m_tools->get_trk_map_node_name());
  auto* sic = findNode::getClass<TrackSeedContainer>(topNode, "SiliconTrackSeedContainer");
  if (!vmap || !tmap)
  {
    return;
  }
  SvtxVertex* sv = vmap->get(pvIn.Id());
  if (!sv)
  {
    return;
  }
  std::set<unsigned int> dids;
  for (int i = 0; i < nd; ++i)
  {
    dids.insert(daughters[i].Id());
    m_dr_inPV[i] = (sv->find_track(daughters[i].Id()) != sv->end_tracks()) ? 1 : 0;
  }
  std::vector<PVfixTrack> all, noD;
  int nNoD = 0;
  for (auto it = sv->begin_tracks(); it != sv->end_tracks(); ++it)
  {
    const bool isD = dids.count(*it) > 0;
    nNoD += isD ? 0 : 1;
    SvtxTrack* t = tmap->get(*it);
    if (!t || t->get_quality() > 300 || t->get_pt() < 0.1)
    {
      continue;
    }
    TrackSeed* ss = t->get_silicon_seed();
    if (!ss && sic)
    {
      ss = KFParticle_Tools::findSeed(sic, t->get_silicon_seed_id());
    }
    int nm = 0;
    if (ss)
    {
      for (auto k = ss->begin_cluster_keys(); k != ss->end_cluster_keys(); ++k)
      {
        nm += (TrkrDefs::getTrkrId(*k) == TrkrDefs::mvtxId) ? 1 : 0;
      }
    }
    if (nm < 3)
    {
      continue;
    }
    PVfixTrack r{};
    r.id = *it;
    r.a[0] = t->get_x();
    r.a[1] = t->get_y();
    r.a[2] = t->get_z();
    double p = t->get_p();
    r.b[0] = t->get_px() / p;
    r.b[1] = t->get_py() / p;
    r.b[2] = t->get_pz() / p;
    for (int a = 0; a < 3; ++a)
    {
      for (int b = 0; b < 3; ++b)
      {
        r.cov[a][b] = t->get_error(a, b);
      }
    }
    all.push_back(r);
    if (!isD)
    {
      noD.push_back(r);
    }
  }
  m_nTracksOfVertex_noDaughters = nNoD;

  // PV-dependent quantities of the candidate w.r.t. a given PV
  auto makePV = [&](const double* ps, const double* cv) {
    KFParticle v;
    float par[6] = {(float) ps[0], (float) ps[1], (float) ps[2], 0, 0, 0};
    float c21[21] = {0};
    for (int k = 0; k < 6; ++k)
    {
      c21[k] = cv[k];
    }
    v.Create(par, c21, 0, -1);
    v.NDF() = 0;
    v.Chi2() = 0;
    v.SetId(pvIn.Id());
    return v;
  };
  auto fillQ = [&](const KFParticle& v, float* mq, float* dca, float* dcaxy, float* sig, float* sigxy) {
    KFParticle m = mother;
    mq[0] = m_tools->eventDIRA(m, v);
    mq[1] = m_tools->eventDIRA(m, v, false);
    mq[2] = m_tools->flightDistanceChi2(m, v);
    mq[3] = m.GetDistanceFromVertex(v);
    mq[4] = m.GetDistanceFromVertexXY(v);
    KFParticle mp = mother;
    mp.SetProductionVertex(v);
    mp.GetDecayLength(mq[5], mq[6]);
    mp.GetDecayLengthXY(mq[7], mq[8]);
    for (int i = 0; i < nd; ++i)
    {
      dca[i] = daughters[i].GetDistanceFromVertex(v);
      dcaxy[i] = daughters[i].GetDistanceFromVertexXY(v);
      sig[i] = daughters[i].GetDeviationFromVertex(v);
      sigxy[i] = daughters[i].GetDeviationFromVertexXY(v);
    }
  };

  double pos[3], cov[6];
  m_pvra_status = pvfixSimpleVertex(all, pos, cov, m_pvra_npair, m_pvra_ntrk);
  if (m_pvra_status == 1)
  {
    m_pvra_x = pos[0];
    m_pvra_y = pos[1];
    m_pvra_z = pos[2];
    float mq[9];
    fillQ(makePV(pos, cov), mq, m_da_dca, m_da_dca_xy, m_da_dca_sig, m_da_dca_sig_xy);
    m_ma_dira = mq[0];
    m_ma_dira_xy = mq[1];
    m_ma_fdchi2 = mq[2];
    m_ma_pvdca = mq[3];
    m_ma_pvdca_xy = mq[4];
    m_ma_L = mq[5];
    m_ma_Lerr = mq[6];
    m_ma_Lxy = mq[7];
    m_ma_Lxyerr = mq[8];
  }

  KFParticle pv;
  const char* closure = std::getenv("KFP_PVREFIT_CLOSURE");
  if (closure && std::string(closure) == "1")
  {
    pv = pvIn;  // code-path closure: the *_PVrefit quantities must reproduce the original ones
    m_pvr_status = 2;
    pos[0] = pvIn.GetX();
    pos[1] = pvIn.GetY();
    pos[2] = pvIn.GetZ();
    for (int k = 0; k < 6; ++k)
    {
      cov[k] = pvIn.GetCovariance(k);
    }
  }
  else
  {
    m_pvr_status = pvfixSimpleVertex(noD, pos, cov, m_pvr_npair, m_pvr_ntrk);
    if (m_pvr_status != 1)
    {
      return;
    }
    pv = makePV(pos, cov);
  }
  m_pvr_x = pos[0];
  m_pvr_y = pos[1];
  m_pvr_z = pos[2];
  for (int k = 0; k < 6; ++k)
  {
    m_pvr_cov[k] = cov[k];
  }

  float mq[9];
  fillQ(pv, mq, m_dr_dca, m_dr_dca_xy, m_dr_dca_sig, m_dr_dca_sig_xy);
  m_mr_dira = mq[0];
  m_mr_dira_xy = mq[1];
  m_mr_fdchi2 = mq[2];
  m_mr_pvdca = mq[3];
  m_mr_pvdca_xy = mq[4];
  m_mr_L = mq[5];
  m_mr_Lerr = mq[6];
  m_mr_Lxy = mq[7];
  m_mr_Lxyerr = mq[8];
}
'''.replace('\n', nl)
wr('KFParticle_nTuple.cc', c)
print('fix2 applied')
