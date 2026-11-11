#include "TrackQATree.h"

#include <ffaobjects/EventHeader.h>
#include <ffarawobjects/Gl1Packet.h>

#include <globalvertex/SvtxVertex.h>
#include <globalvertex/SvtxVertexMap.h>

#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/SvtxTrack_v5.h>
#include <trackbase_historic/TrackAnalysisUtils.h>
#include <trackbase_historic/TrackSeed.h>
#include <trackbase_historic/TrackSeedContainer.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/getClass.h>

#include <TFile.h>
#include <TTree.h>

#include <climits>
#include <cmath>
#include <iostream>
#include <set>

TrackQATree::TrackQATree(const std::string &name, const std::string &outfile)
  : SubsysReco(name)
  , m_outfileName(outfile)
{
}

int TrackQATree::Init(PHCompositeNode * /*topNode*/)
{
  m_file = TFile::Open(m_outfileName.c_str(), "RECREATE");
  if (!m_file || m_file->IsZombie())
  {
    std::cout << Name() << ": cannot open " << m_outfileName << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  m_file->SetCompressionSettings(505);  // ZSTD level 5

  m_evtTree = new TTree("event", "TrackQA: one entry per event");
  m_evtTree->Branch("run", &m_run, "run/I");
  m_evtTree->Branch("event", &m_evt, "event/I");
  m_evtTree->Branch("GL1_BCO", &m_gl1_bco, "GL1_BCO/L");
  m_evtTree->Branch("GL1_scaled", &m_gl1_scaled, "GL1_scaled/l");
  m_evtTree->Branch("GL1_live", &m_gl1_live, "GL1_live/l");
  m_evtTree->Branch("GL1_bunch", &m_gl1_bunch, "GL1_bunch/I");
  m_evtTree->Branch("nTrk", &m_ntrk, "nTrk/I");
  m_evtTree->Branch("nTrkSil", &m_ntrk_sil, "nTrkSil/I");
  m_evtTree->Branch("nTrkINTT", &m_ntrk_intt, "nTrkINTT/I");
  m_evtTree->Branch("nTrkMVTXonly", &m_ntrk_mvtxonly, "nTrkMVTXonly/I");
  m_evtTree->Branch("nTrkTPConly", &m_ntrk_tpconly, "nTrkTPConly/I");
  m_evtTree->Branch("nTrkX0", &m_ntrk_x0, "nTrkX0/I");
  m_evtTree->Branch("nTrkSilX0", &m_ntrk_silx0, "nTrkSilX0/I");
  m_evtTree->Branch("nXingTrk", &m_nxing_trk, "nXingTrk/I");
  m_evtTree->Branch("nVtx", &m_nvtx, "nVtx/I");
  m_evtTree->Branch("nVtxX0", &m_nvtx_x0, "nVtxX0/I");
  m_evtTree->Branch("vx", &m_vx, "vx/F");
  m_evtTree->Branch("vy", &m_vy, "vy/F");
  m_evtTree->Branch("vz", &m_vz, "vz/F");
  m_evtTree->Branch("vntrk", &m_vntrk, "vntrk/I");
  m_evtTree->Branch("vchi2ndf", &m_vchi2ndf, "vchi2ndf/F");

  m_trkTree = new TTree("track", "TrackQA: one entry per track");
  m_trkTree->Branch("run", &m_run, "run/I");
  m_trkTree->Branch("event", &m_evt, "event/I");
  m_trkTree->Branch("ievt", &m_t_ievt, "ievt/L");
  m_trkTree->Branch("id", &m_t_id, "id/I");
  m_trkTree->Branch("pt", &m_pt, "pt/F");
  m_trkTree->Branch("eta", &m_eta, "eta/F");
  m_trkTree->Branch("phi", &m_phi, "phi/F");
  m_trkTree->Branch("charge", &m_charge, "charge/S");
  m_trkTree->Branch("pcax", &m_pcax, "pcax/F");
  m_trkTree->Branch("pcay", &m_pcay, "pcay/F");
  m_trkTree->Branch("pcaz", &m_pcaz, "pcaz/F");
  m_trkTree->Branch("chi2", &m_chi2, "chi2/F");
  m_trkTree->Branch("ndf", &m_ndf, "ndf/S");
  m_trkTree->Branch("quality", &m_quality, "quality/F");
  m_trkTree->Branch("crossing", &m_crossing, "crossing/S");
  m_trkTree->Branch("nmvtx", &m_nmvtx, "nmvtx/b");
  m_trkTree->Branch("nintt", &m_nintt, "nintt/b");
  m_trkTree->Branch("ntpc", &m_ntpc, "ntpc/b");
  m_trkTree->Branch("ntpot", &m_ntpot, "ntpot/b");
  m_trkTree->Branch("hasSiSeed", &m_hasSiSeed, "hasSiSeed/O");
  m_trkTree->Branch("hasTpcSeed", &m_hasTpcSeed, "hasTpcSeed/O");
  m_trkTree->Branch("vtxMatch", &m_vtxMatch, "vtxMatch/S");
  m_trkTree->Branch("vtxId", &m_vtxId, "vtxId/I");
  m_trkTree->Branch("tvx", &m_tvx, "tvx/F");
  m_trkTree->Branch("tvy", &m_tvy, "tvy/F");
  m_trkTree->Branch("tvz", &m_tvz, "tvz/F");
  m_trkTree->Branch("tvntrk", &m_tvntrk, "tvntrk/S");
  m_trkTree->Branch("dcaxy", &m_dcaxy, "dcaxy/F");
  m_trkTree->Branch("dcaxyErr", &m_dcaxyErr, "dcaxyErr/F");
  m_trkTree->Branch("dcaz", &m_dcaz, "dcaz/F");
  m_trkTree->Branch("dcazErr", &m_dcazErr, "dcazErr/F");
  m_trkTree->Branch("tpcSeedPt", &m_tpcSeedPt, "tpcSeedPt/F");
  m_trkTree->Branch("tpcSeedEta", &m_tpcSeedEta, "tpcSeedEta/F");
  m_trkTree->Branch("tpcSeedPhi", &m_tpcSeedPhi, "tpcSeedPhi/F");
  return Fun4AllReturnCodes::EVENT_OK;
}

int TrackQATree::process_event(PHCompositeNode *topNode)
{
  m_run = m_evt = -1;
  m_gl1_bco = -1;
  m_gl1_scaled = m_gl1_live = 0;
  m_gl1_bunch = -1;
  m_ntrk = m_ntrk_sil = m_ntrk_intt = m_ntrk_mvtxonly = m_ntrk_tpconly = m_ntrk_x0 = m_ntrk_silx0 = 0;
  m_nvtx = m_nvtx_x0 = m_nxing_trk = 0;
  m_vx = m_vy = m_vz = m_vchi2ndf = -999;
  m_vntrk = 0;

  if (auto *evth = findNode::getClass<EventHeader>(topNode, "EventHeader"))
  {
    m_run = evth->get_RunNumber();
    m_evt = evth->get_EvtSequence();
  }
  auto *gl1 = findNode::getClass<Gl1Packet>(topNode, "GL1RAWHIT");
  if (!gl1)
  {
    gl1 = findNode::getClass<Gl1Packet>(topNode, "GL1Packet");
  }
  if (gl1)
  {
    m_gl1_bco = gl1->lValue(0, "BCO");
    m_gl1_scaled = gl1->getScaledVector();
    m_gl1_live = gl1->getLiveVector();
    m_gl1_bunch = gl1->getBunchNumber();
  }

  // ---- vertices
  auto *vtxmap = findNode::getClass<SvtxVertexMap>(topNode, m_vertexMapName);
  if (vtxmap)
  {
    for (const auto &[id, vtx] : *vtxmap)
    {
      ++m_nvtx;
      if (vtx->get_beam_crossing() != 0)
      {
        continue;
      }
      ++m_nvtx_x0;
      if (static_cast<int>(vtx->size_tracks()) > m_vntrk)
      {
        m_vntrk = vtx->size_tracks();
        m_vx = vtx->get_x();
        m_vy = vtx->get_y();
        m_vz = vtx->get_z();
        m_vchi2ndf = vtx->get_ndof() > 0 ? vtx->get_chisq() / vtx->get_ndof() : -999;
      }
    }
  }

  // ---- tracks
  auto *tpcSeeds = findNode::getClass<TrackSeedContainer>(topNode, m_tpcSeedName);
  auto *trkmap = findNode::getClass<SvtxTrackMap>(topNode, m_trackMapName);
  std::set<int> xings;
  if (trkmap)
  {
    for (const auto &[id, trk] : *trkmap)
    {
      ++m_ntrk;
      int nmvtx = 0, nintt = 0, ntpc = 0, ntpot = 0;
      if (auto *v5 = dynamic_cast<SvtxTrack_v5 *>(trk))
      {
        nmvtx = v5->get_nmvtx_clusters();
        nintt = v5->get_nintt_clusters();
        ntpc = v5->get_ntpc_clusters();
        ntpot = v5->get_ntpot_clusters();
      }
      else
      {
        auto c = TrackAnalysisUtils::get_cluster_counts(trk);
        nmvtx = std::get<0>(c);
        nintt = std::get<1>(c);
        ntpc = std::get<2>(c);
        ntpot = std::get<3>(c);
      }
      const int xing = trk->get_crossing();
      const bool sil = (nmvtx + nintt) > 0;
      if (sil)
      {
        ++m_ntrk_sil;
        xings.insert(xing);
        if (xing == 0) ++m_ntrk_silx0;
      }
      if (nintt > 0) ++m_ntrk_intt;
      else if (nmvtx > 0) ++m_ntrk_mvtxonly;
      else ++m_ntrk_tpconly;
      if (xing == 0) ++m_ntrk_x0;

      if (trk->get_pt() < m_minPt)
      {
        continue;
      }
      m_t_ievt = m_ievt;
      m_t_id = id;
      m_pt = trk->get_pt();
      m_eta = trk->get_eta();
      m_phi = trk->get_phi();
      m_charge = trk->get_charge();
      m_pcax = trk->get_x();
      m_pcay = trk->get_y();
      m_pcaz = trk->get_z();
      m_chi2 = trk->get_chisq();
      m_ndf = trk->get_ndf();
      m_quality = trk->get_quality();
      m_crossing = xing;
      m_nmvtx = nmvtx;
      m_nintt = nintt;
      m_ntpc = ntpc;
      m_ntpot = ntpot;
      // seed pointers or seed ids (whichever the stored track version provides)
      m_hasSiSeed = trk->get_silicon_seed() != nullptr || trk->get_silicon_seed_id() != USHRT_MAX;
      m_hasTpcSeed = trk->get_tpc_seed() != nullptr || trk->get_tpc_seed_id() != USHRT_MAX;

      // TPC seed kinematics (pointer if set, else id into the TPC seed container)
      m_tpcSeedPt = m_tpcSeedEta = m_tpcSeedPhi = -999;
      {
        TrackSeed *ts = trk->get_tpc_seed();
        if (!ts && tpcSeeds && trk->get_tpc_seed_id() != USHRT_MAX) ts = tpcSeeds->get(trk->get_tpc_seed_id());
        if (ts)
        {
          m_tpcSeedPt = ts->get_pt();
          m_tpcSeedEta = ts->get_eta();
          m_tpcSeedPhi = ts->get_phi();
        }
      }

      // vertex association (see header)
      SvtxVertex *vtx = nullptr;
      m_vtxMatch = 0;
      m_vtxId = -1;
      if (vtxmap)
      {
        const unsigned int vid = trk->get_vertex_id();
        if (vid != UINT_MAX && vtxmap->find(vid) != vtxmap->end())
        {
          vtx = vtxmap->get(vid);
          m_vtxId = vid;
          m_vtxMatch = 2;
          for (auto it = vtx->begin_tracks(); it != vtx->end_tracks(); ++it)
          {
            if (*it == id)
            {
              m_vtxMatch = 1;  // track used in the vertex fit
              break;
            }
          }
        }
        else
        {
          float best = 1e9;
          for (const auto &[vid2, v] : *vtxmap)
          {
            if (static_cast<int>(v->get_beam_crossing()) != xing) continue;
            const float dz = std::fabs(v->get_z() - m_pcaz);
            if (dz < best)
            {
              best = dz;
              vtx = v;
              m_vtxId = vid2;
            }
          }
          if (vtx) m_vtxMatch = 3;
        }
      }
      ++m_nMatch[m_vtxMatch];
      if (vtx)
      {
        m_tvx = vtx->get_x();
        m_tvy = vtx->get_y();
        m_tvz = vtx->get_z();
        m_tvntrk = vtx->size_tracks();
        Acts::Vector3 pv(vtx->get_x(), vtx->get_y(), vtx->get_z());
        auto dca = TrackAnalysisUtils::get_dca(trk, pv);
        m_dcaxy = dca.first.first;
        m_dcaxyErr = dca.first.second;
        m_dcaz = dca.second.first;
        m_dcazErr = dca.second.second;
      }
      else
      {
        m_tvx = m_tvy = m_tvz = -999;
        m_tvntrk = 0;
        m_dcaxy = m_dcaxyErr = m_dcaz = m_dcazErr = -999;
      }
      m_trkTree->Fill();
      ++m_nTrkWritten;
    }
  }
  m_nxing_trk = xings.size();

  m_evtTree->Fill();
  ++m_ievt;
  return Fun4AllReturnCodes::EVENT_OK;
}

int TrackQATree::End(PHCompositeNode * /*topNode*/)
{
  if (m_file)
  {
    m_file->cd();
    m_evtTree->Write();
    m_trkTree->Write();
    m_file->Close();
    std::cout << Name() << ": events " << m_ievt << " tracks written " << m_nTrkWritten
              << " vtxMatch(0,1,2,3) " << m_nMatch[0] << " " << m_nMatch[1] << " " << m_nMatch[2] << " " << m_nMatch[3] << std::endl;
    std::cout << Name() << ": wrote " << m_outfileName << std::endl;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}
