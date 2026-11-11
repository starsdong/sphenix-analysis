#!/usr/bin/env python3
"""d0-pvfix fix 4: selection with the refitted PVs at fill time.
New branches (always filled):
  pvfix_passOriginal : the stock in-job selection, reproduced: both daughters have minimum 3D DCA over the PVs of their
                       crossing (KFParticle_Tools::calcMinPV_DCA, as isGoodTrack) >= dcaCut, and DIRA (original PV) >= diraCut
  pvfix_passRefit    : both |DCA| >= dcaCut and DIRA >= diraCut w.r.t. the pair refit w/o daughters (primary_vertex_refit)
  pvfix_passKF       : the same w.r.t. the KFParticle fit w/o daughters (primary_vertex_kf)
With KFP_PVFIX_POSTCUT="<dcaCut>,<diraCut>" in the environment (cm, cosine), a candidate is written only if
pvfix_passOriginal || pvfix_passRefit || pvfix_passKF; the in-job cuts can then be loosened in the macro so that
candidates biased by the original PV are recovered. Without the variable nothing is rejected and the flags use
dcaCut = diraCut = -inf. Requires fixes 2 and 3. Run in offline/packages/KFParticle_sPHENIX."""

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
h = sub(h, '''  PVfixV m_v_kf, m_v_kfAll;
''', '''  PVfixV m_v_kf, m_v_kfAll;

  // d0-pvfix fix 4: selection with the refitted PVs at fill time
  int m_pvfix_passOriginal{1};
  int m_pvfix_passRefit{0};
  int m_pvfix_passKF{0};
  bool m_pvfix_keep{true};
  bool m_pvfix_post{false};
''', nl)
wr('KFParticle_nTuple.h', h)

c, nl = rd('KFParticle_nTuple.cc')
c = sub(c, '''    bookV("primary_vertex_kf", m_v_kf);
''', '''    bookI("pvfix_passOriginal", &m_pvfix_passOriginal);
    bookI("pvfix_passRefit", &m_pvfix_passRefit);
    bookI("pvfix_passKF", &m_pvfix_passKF);
    bookV("primary_vertex_kf", m_v_kf);
''', nl)
# at the start of pvRefit: defaults
c = sub(c, '''  m_nTracksOfVertex_noDaughters = -1;
  for (auto* q : {&m_q_kf, &m_q_kfAll})''', '''  m_nTracksOfVertex_noDaughters = -1;
  m_pvfix_passOriginal = 1;
  m_pvfix_passRefit = m_pvfix_passKF = 0;
  m_pvfix_keep = true;
  for (auto* q : {&m_q_kf, &m_q_kfAll})''', nl)
# the decision: computed at the very end of pvRefit and also when it returns early -> use a scope guard
c = sub(c, '''  auto* vmap = findNode::getClass<SvtxVertexMap>(topNode, m_tools->get_vtx_map_node_name());
  auto* tmap = findNode::getClass<SvtxTrackMap>(topNode, m_tools->get_trk_map_node_name());''', '''  // d0-pvfix fix 4: post-selection cuts from the environment (parsed once)
  static float s_dcaCut = -std::numeric_limits<float>::max(), s_diraCut = -std::numeric_limits<float>::max();
  static bool s_post = [] {
    const char* e = std::getenv("KFP_PVFIX_POSTCUT");
    if (!e)
    {
      return false;
    }
    if (std::sscanf(e, "%f,%f", &s_dcaCut, &s_diraCut) != 2)
    {
      std::cout << "KFP_PVFIX_POSTCUT: cannot parse '" << e << "', post-selection off" << std::endl;
      s_dcaCut = s_diraCut = -std::numeric_limits<float>::max();
      return false;
    }
    std::cout << "KFP_PVFIX_POSTCUT: keep candidates passing |DCA| >= " << s_dcaCut << " cm and DIRA >= " << s_diraCut
              << " with the original PV (min over PVs, as in-job) OR the pair refit OR the KF refit w/o daughters" << std::endl;
    return true;
  }();
  m_pvfix_post = s_post;
  {
    // original selection as in KFParticle_Tools::isGoodTrack (minimum DCA over the PVs of the track's crossing) + DIRA
    std::vector<KFParticle> allPVs = m_tools->makeAllPrimaryVertices(topNode, m_tools->get_vtx_map_node_name());
    bool pass = m_calculated_mother_dira >= s_diraCut;
    for (int i = 0; i < nd; ++i)
    {
      float dmin = 0, smin = 0;
      m_tools->calcMinPV_DCA(daughters[i], allPVs, dmin, smin);
      pass = pass && (dmin >= s_dcaCut);
    }
    m_pvfix_passOriginal = pass ? 1 : 0;
  }

  auto* vmap = findNode::getClass<SvtxVertexMap>(topNode, m_tools->get_vtx_map_node_name());
  auto* tmap = findNode::getClass<SvtxTrackMap>(topNode, m_tools->get_trk_map_node_name());''', nl)
# KF pass flag, after doKF(noD, ...)
c = sub(c, '''  doKF(noD, m_v_kf, m_q_kf);
''', '''  doKF(noD, m_v_kf, m_q_kf);
  if (m_v_kf.status == 1)
  {
    bool pass = m_q_kf.m[0] >= s_diraCut;
    for (int i = 0; i < nd; ++i)
    {
      pass = pass && (m_q_kf.dca[i] >= s_dcaCut);
    }
    m_pvfix_passKF = pass ? 1 : 0;
  }
''', nl)
# refit pass flag, at the end of pvRefit (after the refit quantities are filled)
c = sub(c, '''  float mq[9];
  fillQ(pv, mq, m_dr_dca, m_dr_dca_xy, m_dr_dca_sig, m_dr_dca_sig_xy);
''', '''  float mq[9];
  fillQ(pv, mq, m_dr_dca, m_dr_dca_xy, m_dr_dca_sig, m_dr_dca_sig_xy);
  {
    bool pass = mq[0] >= s_diraCut;
    for (int i = 0; i < nd; ++i)
    {
      pass = pass && (m_dr_dca[i] >= s_dcaCut);
    }
    m_pvfix_passRefit = (m_pvr_status == 1 && pass) ? 1 : 0;
  }
''', nl)
# keep decision right after the refit call in fillBranch
c = sub(c, '''    pvRefit(topNode, motherForRefit, vertex_fillbranch, daughters);  // d0-pvfix fix 2
''', '''    pvRefit(topNode, motherForRefit, vertex_fillbranch, daughters);  // d0-pvfix fix 2
    m_pvfix_keep = !m_pvfix_post || m_pvfix_passOriginal || m_pvfix_passRefit || m_pvfix_passKF;  // d0-pvfix fix 4
''', nl)
# write only kept candidates
c = sub(c, '''  if (fillConditionMet())
  {
    m_tree->Fill();
  }
''', '''  if (fillConditionMet() && m_pvfix_keep)  // d0-pvfix fix 4
  {
    m_tree->Fill();
  }
''', nl)
c = sub(c, '''#include <array>
#include <limits>
#include <memory>
''', '''#include <array>
#include <cstdio>
#include <iostream>
#include <limits>
#include <memory>
''', nl)
wr('KFParticle_nTuple.cc', c)
print('fix4 applied')
