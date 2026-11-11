// EventInfoTree: per-event companion tree for the KFParticle ntuples.
//
// One entry per processed event (every event, also those without candidates), so
//  - yields can be normalised to the number of analysed events,
//  - the crossing of the PV picked by KFParticle can be looked up via
//    (runNumber, eventNumber, primary_vertex_ID) -> vtx_crossing,
//  - the track crossing population of the event (silicon / INTT / TPC-only) is known.
#ifndef EVENTINFOTREE_H
#define EVENTINFOTREE_H

#include <fun4all/SubsysReco.h>

#include <cstdint>
#include <string>
#include <vector>

class PHCompositeNode;
class TFile;
class TTree;

class EventInfoTree : public SubsysReco
{
 public:
  explicit EventInfoTree(const std::string &name = "EventInfoTree", const std::string &outfile = "eventinfo.root");
  ~EventInfoTree() override = default;

  int Init(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;
  int End(PHCompositeNode *topNode) override;

  void setTrackMapName(const std::string &n) { m_trackMapName = n; }
  void setVertexMapName(const std::string &n) { m_vertexMapName = n; }

 private:
  std::string m_outfileName;
  std::string m_trackMapName{"SvtxTrackMap"};
  std::string m_vertexMapName{"SvtxVertexMap"};
  TFile *m_file{nullptr};
  TTree *m_tree{nullptr};

  // event level
  int m_run{-1};
  int m_evt{-1};
  int64_t m_gl1_bco{-1};
  uint64_t m_gl1_scaled{0};
  uint64_t m_gl1_live{0};
  int m_gl1_bunch{-1};
  // track counts
  int m_ntrk{0};          // all tracks in SvtxTrackMap
  int m_ntrk_intt{0};     // >=1 INTT cluster (crossing from INTT timing)
  int m_ntrk_mvtxonly{0}; // MVTX but no INTT clusters
  int m_ntrk_tpconly{0};  // no silicon clusters (crossing defaults to 0)
  int m_ntrk_x0{0};       // tracks with crossing == 0
  // vertices
  int m_nvtx{0};
  std::vector<int> m_vtx_id;
  std::vector<int> m_vtx_crossing;
  std::vector<float> m_vtx_x, m_vtx_y, m_vtx_z;
  std::vector<float> m_vtx_chi2;
  std::vector<int> m_vtx_ndof;
  std::vector<int> m_vtx_ntrk;
};

#endif
