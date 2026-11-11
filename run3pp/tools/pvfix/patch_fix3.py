#!/usr/bin/env python3
"""d0-pvfix fix 3: KFParticle vertex fit (KFVertex::ConstructPrimaryVertex, chi2 cut 3.5) of the tracks attached to the
candidate's SvtxVertex, with all tracks (primary_vertex_kfAll_*, *_PVkfAll) and without the candidate's daughters
(primary_vertex_kf_*, *_PVkf). Same track list and quality cuts as the fix-2 refits; the tracks enter with their full
DST state (position, momentum, 6x6 covariance), like KFParticle_Tools::makeParticle. Requires fix 2. Run in
offline/packages/KFParticle_sPHENIX."""

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
h = sub(h, '''  float m_da_dca_sig_xy[kMaxRefitTracks]{0};
''', '''  float m_da_dca_sig_xy[kMaxRefitTracks]{0};

  // d0-pvfix fix 3: KFParticle vertex fit with all attached tracks (kfAll) and without the daughters (kf)
  struct PVfixQ
  {
    float m[9]{0};  // DIRA, DIRA_xy, FDchi2, PV_DCA, PV_DCA_xy, decayLength, decayLengthErr, decayLength_xy, decayLengthErr_xy
    float dca[kMaxRefitTracks]{0};
    float dcaxy[kMaxRefitTracks]{0};
    float sig[kMaxRefitTracks]{0};
    float sigxy[kMaxRefitTracks]{0};
  };
  struct PVfixV
  {
    int status{-1};
    int ntrk{-1};
    int ndf{-1};
    float x{0};
    float y{0};
    float z{0};
    float chi2{0};
  };
  PVfixQ m_q_kf, m_q_kfAll;
  PVfixV m_v_kf, m_v_kfAll;
''', nl)
wr('KFParticle_nTuple.h', h)

c, nl = rd('KFParticle_nTuple.cc')
# booking
c = sub(c, '''    for (int i = 0; i < std::min(m_num_tracks_nTuple, kMaxRefitTracks); ++i)
    {
      std::string dn = "track_" + std::to_string(i + 1);
      bookI(dn + "_inPV", &m_dr_inPV[i]);
''', '''    // d0-pvfix fix 3
    auto bookV = [&](const std::string& n, PVfixV& v) {
      bookF(n + "_x", &v.x);
      bookF(n + "_y", &v.y);
      bookF(n + "_z", &v.z);
      bookF(n + "_chi2", &v.chi2);
      bookI(n + "_nDoF", &v.ndf);
      bookI(n + "_status", &v.status);
      bookI(n + "_nTracks", &v.ntrk);
    };
    auto bookQ = [&](const std::string& sfx, PVfixQ& q) {
      const char* mn[9] = {"_DIRA", "_DIRA_xy", "_FDchi2", "_PV_DCA", "_PV_DCA_xy", "_decayLength", "_decayLengthErr", "_decayLength_xy", "_decayLengthErr_xy"};
      for (int k = 0; k < 9; ++k)
      {
        bookF(mother_name + mn[k] + sfx, &q.m[k]);
      }
      for (int i = 0; i < std::min(m_num_tracks_nTuple, kMaxRefitTracks); ++i)
      {
        std::string dn = "track_" + std::to_string(i + 1);
        bookF(dn + "_PV_DCA" + sfx, &q.dca[i]);
        bookF(dn + "_PV_DCA_xy" + sfx, &q.dcaxy[i]);
        bookF(dn + "_PV_DCA_sig" + sfx, &q.sig[i]);
        bookF(dn + "_PV_DCA_sig_xy" + sfx, &q.sigxy[i]);
      }
    };
    bookV("primary_vertex_kf", m_v_kf);
    bookV("primary_vertex_kfAll", m_v_kfAll);
    bookQ("_PVkf", m_q_kf);
    bookQ("_PVkfAll", m_q_kfAll);
    for (int i = 0; i < std::min(m_num_tracks_nTuple, kMaxRefitTracks); ++i)
    {
      std::string dn = "track_" + std::to_string(i + 1);
      bookI(dn + "_inPV", &m_dr_inPV[i]);
''', nl)

# track pointer in PVfixTrack
c = sub(c, '''  struct PVfixTrack
  {
    unsigned int id;
''', '''  struct PVfixTrack
  {
    SvtxTrack* trk;
    unsigned int id;
''', nl)
c = sub(c, '''    PVfixTrack r{};
    r.id = *it;
''', '''    PVfixTrack r{};
    r.trk = t;
    r.id = *it;
''', nl)

# KF fit helper before pvRefit
c = sub(c, '''}  // namespace

void KFParticle_nTuple::pvRefit(''', '''  // d0-pvfix fix 3: KFParticle primary-vertex fit of a set of tracks; returns 1 if >= 2 tracks remain
  int pvfixKFVertex(const std::vector<PVfixTrack>& t, KFParticle& out, int& nUsed)
  {
    nUsed = 0;
    if (t.size() < 2)
    {
      return 0;
    }
    std::vector<KFParticle> parts(t.size());
    std::vector<const KFParticle*> ptr(t.size());
    std::unique_ptr<Bool_t[]> flag(new Bool_t[t.size()]);
    for (size_t i = 0; i < t.size(); ++i)
    {
      const SvtxTrack* s = t[i].trk;
      float par[6] = {s->get_x(), s->get_y(), s->get_z(), s->get_px(), s->get_py(), s->get_pz()};
      float cov[21];
      int k = 0;
      for (int a = 0; a < 6; ++a)
      {
        for (int b = 0; b <= a; ++b)
        {
          cov[k++] = s->get_error(a, b);
        }
      }
      parts[i].Create(par, cov, (Int_t) s->get_charge(), -1);
      parts[i].NDF() = s->get_ndf();
      parts[i].Chi2() = s->get_chisq();
      parts[i].SetId(s->get_id());
      ptr[i] = &parts[i];
      flag[i] = true;
    }
    KFVertex v;
    v.ConstructPrimaryVertex(ptr.data(), (int) t.size(), flag.get(), 3.5);
    for (size_t i = 0; i < t.size(); ++i)
    {
      nUsed += flag[i] ? 1 : 0;
    }
    if (nUsed < 2 || !std::isfinite(v.GetX()) || !std::isfinite(v.GetZ()))
    {
      return 0;
    }
    out = v;
    return 1;
  }
}  // namespace

void KFParticle_nTuple::pvRefit(''', nl)

# resets
c = sub(c, '''  m_nTracksOfVertex_noDaughters = -1;
''', '''  m_nTracksOfVertex_noDaughters = -1;
  for (auto* q : {&m_q_kf, &m_q_kfAll})
  {
    for (float& x : q->m)
    {
      x = nan;
    }
    for (int i = 0; i < kMaxRefitTracks; ++i)
    {
      q->dca[i] = q->dcaxy[i] = q->sig[i] = q->sigxy[i] = nan;
    }
  }
  for (auto* v : {&m_v_kf, &m_v_kfAll})
  {
    v->status = v->ntrk = v->ndf = -1;
    v->x = v->y = v->z = v->chi2 = nan;
  }
''', nl)

# compute the KF versions right after the all-tracks pair refit block (fillQ and makePV exist there)
c = sub(c, '''  KFParticle pv;
  const char* closure = std::getenv("KFP_PVREFIT_CLOSURE");''', '''  // d0-pvfix fix 3: KFParticle vertex fits
  auto doKF = [&](const std::vector<PVfixTrack>& trks, PVfixV& v, PVfixQ& q) {
    KFParticle kv;
    v.status = pvfixKFVertex(trks, kv, v.ntrk);
    if (v.status != 1)
    {
      return;
    }
    v.x = kv.GetX();
    v.y = kv.GetY();
    v.z = kv.GetZ();
    v.chi2 = kv.GetChi2();
    v.ndf = kv.GetNDF();
    kv.SetId(pvIn.Id());
    fillQ(kv, q.m, q.dca, q.dcaxy, q.sig, q.sigxy);
  };
  doKF(all, m_v_kfAll, m_q_kfAll);
  doKF(noD, m_v_kf, m_q_kf);

  KFParticle pv;
  const char* closure = std::getenv("KFP_PVREFIT_CLOSURE");''', nl)
c = sub(c, '''#include <array>
#include <limits>
#include <set>
''', '''#include <array>
#include <limits>
#include <memory>
#include <set>
''', nl)
wr('KFParticle_nTuple.cc', c)
print('fix3 applied')
