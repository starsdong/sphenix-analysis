#include "EventInfoTree.h"

#include <ffaobjects/EventHeader.h>
#include <ffarawobjects/Gl1Packet.h>

#include <globalvertex/SvtxVertex.h>
#include <globalvertex/SvtxVertexMap.h>

#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/SvtxTrack_v5.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/getClass.h>

#include <TFile.h>
#include <TTree.h>

#include <iostream>

EventInfoTree::EventInfoTree(const std::string &name, const std::string &outfile)
  : SubsysReco(name)
  , m_outfileName(outfile)
{
}

int EventInfoTree::Init(PHCompositeNode * /*topNode*/)
{
  m_file = TFile::Open(m_outfileName.c_str(), "RECREATE");
  if (!m_file || m_file->IsZombie())
  {
    std::cout << Name() << ": cannot open " << m_outfileName << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  m_file->SetCompressionSettings(505);  // ZSTD level 5
  m_tree = new TTree("EventInfo", "per-event info (all processed events)");
  m_tree->Branch("runNumber", &m_run, "runNumber/I");
  m_tree->Branch("eventNumber", &m_evt, "eventNumber/I");
  m_tree->Branch("GL1_BCO", &m_gl1_bco, "GL1_BCO/L");
  m_tree->Branch("GL1_scaled", &m_gl1_scaled, "GL1_scaled/l");
  m_tree->Branch("GL1_live", &m_gl1_live, "GL1_live/l");
  m_tree->Branch("GL1_bunch", &m_gl1_bunch, "GL1_bunch/I");
  m_tree->Branch("nTrk", &m_ntrk, "nTrk/I");
  m_tree->Branch("nTrkINTT", &m_ntrk_intt, "nTrkINTT/I");
  m_tree->Branch("nTrkMVTXonly", &m_ntrk_mvtxonly, "nTrkMVTXonly/I");
  m_tree->Branch("nTrkTPConly", &m_ntrk_tpconly, "nTrkTPConly/I");
  m_tree->Branch("nTrkX0", &m_ntrk_x0, "nTrkX0/I");
  m_tree->Branch("nVtx", &m_nvtx, "nVtx/I");
  m_tree->Branch("vtx_id", &m_vtx_id);
  m_tree->Branch("vtx_crossing", &m_vtx_crossing);
  m_tree->Branch("vtx_x", &m_vtx_x);
  m_tree->Branch("vtx_y", &m_vtx_y);
  m_tree->Branch("vtx_z", &m_vtx_z);
  m_tree->Branch("vtx_chi2", &m_vtx_chi2);
  m_tree->Branch("vtx_ndof", &m_vtx_ndof);
  m_tree->Branch("vtx_ntrk", &m_vtx_ntrk);
  return Fun4AllReturnCodes::EVENT_OK;
}

int EventInfoTree::process_event(PHCompositeNode *topNode)
{
  m_run = m_evt = -1;
  m_gl1_bco = -1;
  m_gl1_scaled = m_gl1_live = 0;
  m_gl1_bunch = -1;
  m_ntrk = m_ntrk_intt = m_ntrk_mvtxonly = m_ntrk_tpconly = m_ntrk_x0 = 0;
  m_nvtx = 0;
  m_vtx_id.clear();
  m_vtx_crossing.clear();
  m_vtx_x.clear();
  m_vtx_y.clear();
  m_vtx_z.clear();
  m_vtx_chi2.clear();
  m_vtx_ndof.clear();
  m_vtx_ntrk.clear();

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

  if (auto *trkmap = findNode::getClass<SvtxTrackMap>(topNode, m_trackMapName))
  {
    for (const auto &[id, trk] : *trkmap)
    {
      ++m_ntrk;
      if (trk->get_crossing() == 0)
      {
        ++m_ntrk_x0;
      }
      int nmvtx = 0, nintt = 0;
      if (auto *v5 = dynamic_cast<SvtxTrack_v5 *>(trk))
      {
        nmvtx = v5->get_nmvtx_clusters();
        nintt = v5->get_nintt_clusters();
      }
      if (nintt > 0)
      {
        ++m_ntrk_intt;
      }
      else if (nmvtx > 0)
      {
        ++m_ntrk_mvtxonly;
      }
      else
      {
        ++m_ntrk_tpconly;
      }
    }
  }

  if (auto *vtxmap = findNode::getClass<SvtxVertexMap>(topNode, m_vertexMapName))
  {
    for (const auto &[id, vtx] : *vtxmap)
    {
      ++m_nvtx;
      m_vtx_id.push_back(id);
      m_vtx_crossing.push_back(vtx->get_beam_crossing());
      m_vtx_x.push_back(vtx->get_x());
      m_vtx_y.push_back(vtx->get_y());
      m_vtx_z.push_back(vtx->get_z());
      m_vtx_chi2.push_back(vtx->get_chisq());
      m_vtx_ndof.push_back(vtx->get_ndof());
      m_vtx_ntrk.push_back(vtx->size_tracks());
    }
  }

  m_tree->Fill();
  return Fun4AllReturnCodes::EVENT_OK;
}

int EventInfoTree::End(PHCompositeNode * /*topNode*/)
{
  if (m_file)
  {
    m_file->cd();
    m_tree->Write();
    m_file->Close();
    std::cout << Name() << ": wrote " << m_outfileName << std::endl;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}
