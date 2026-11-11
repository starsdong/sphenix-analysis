// TrackQATree: trim a DST_TRKR_TRACKS file to a flat event + track ntuple for event/track QA.
//
// Output file (ZSTD) with two trees:
//   "event" : one entry per processed event (all events): run/event, GL1 info, track-type counts,
//             vertex counts, the crossing-0 vertex with most tracks.
//   "track" : one entry per SvtxTrack (all tracks, optional pT threshold): kinematics, fit quality,
//             cluster counts per detector, beam crossing, DCA to the associated vertex.
//             Branches run/event/ievt link a track to its event (ievt = entry number in "event").
//             v2: tpcSeedPt/Eta/Phi = kinematics of the track's TPC seed (TPC-only estimate for every track).
//
// Vertex association for the track DCA (vtxMatch branch):
//   1 : vertex = SvtxTrack::get_vertex_id() and the track is in that vertex's track list (used in the fit)
//   2 : vertex = SvtxTrack::get_vertex_id(), track not used in the vertex fit
//   3 : no vertex id on the track; vertex with the SAME beam crossing closest in z to the track PCA
//   0 : no vertex found -> DCA and vertex branches = -999
// (check on 200 events of 79509-00000, tools/check_vertex_id.C: the vertex id, when set, always
//  points to a vertex of the track's own crossing.)
// TPC-only tracks always carry crossing 0 (no silicon timing), so their DCA is meaningful only if
// they really come from the triggered crossing.
#ifndef TRACKQATREE_H
#define TRACKQATREE_H

#include <fun4all/SubsysReco.h>

#include <cstdint>
#include <string>

class PHCompositeNode;
class TFile;
class TTree;

class TrackQATree : public SubsysReco
{
 public:
  explicit TrackQATree(const std::string &name = "TrackQATree", const std::string &outfile = "trackqa.root");
  ~TrackQATree() override = default;

  int Init(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;
  int End(PHCompositeNode *topNode) override;

  void setTrackMapName(const std::string &n) { m_trackMapName = n; }
  void setVertexMapName(const std::string &n) { m_vertexMapName = n; }
  void setTpcSeedName(const std::string &n) { m_tpcSeedName = n; }
  void setMinPt(float v) { m_minPt = v; }  // tracks below are counted in "event" but not written to "track"

 private:
  std::string m_outfileName;
  std::string m_trackMapName{"SvtxTrackMap"};
  std::string m_vertexMapName{"SvtxVertexMap"};
  std::string m_tpcSeedName{"TpcTrackSeedContainer"};  // absent in trimmed DSTs -> seed branches -999
  float m_minPt{0.};
  TFile *m_file{nullptr};
  TTree *m_evtTree{nullptr};
  TTree *m_trkTree{nullptr};
  long m_ievt{0};
  long m_nTrkWritten{0};
  long m_nMatch[4]{0, 0, 0, 0};

  // ---- event tree
  int m_run{-1};
  int m_evt{-1};
  int64_t m_gl1_bco{-1};
  uint64_t m_gl1_scaled{0};
  uint64_t m_gl1_live{0};
  int m_gl1_bunch{-1};
  int m_ntrk{0};          // all tracks in SvtxTrackMap
  int m_ntrk_sil{0};      // >=1 MVTX or INTT cluster
  int m_ntrk_intt{0};     // >=1 INTT cluster
  int m_ntrk_mvtxonly{0}; // MVTX but no INTT clusters
  int m_ntrk_tpconly{0};  // no silicon clusters (crossing defaults to 0)
  int m_ntrk_x0{0};       // crossing == 0 (incl. TPC-only)
  int m_ntrk_silx0{0};    // silicon tracks with crossing == 0
  int m_nvtx{0};
  int m_nvtx_x0{0};       // vertices with crossing 0
  int m_nxing_trk{0};     // distinct crossings among silicon tracks
  float m_vx{-999}, m_vy{-999}, m_vz{-999};  // crossing-0 vertex with most tracks
  int m_vntrk{0};
  float m_vchi2ndf{-999};

  // ---- track tree
  long m_t_ievt{0};
  int m_t_id{0};
  float m_pt{0}, m_eta{0}, m_phi{0};
  short m_charge{0};
  float m_pcax{0}, m_pcay{0}, m_pcaz{0};
  float m_chi2{0};
  short m_ndf{0};
  float m_quality{0};
  short m_crossing{0};
  unsigned char m_nmvtx{0}, m_nintt{0}, m_ntpc{0}, m_ntpot{0};
  bool m_hasSiSeed{false}, m_hasTpcSeed{false};
  short m_vtxMatch{0};
  int m_vtxId{-1};
  float m_tvx{-999}, m_tvy{-999}, m_tvz{-999};
  short m_tvntrk{0};
  float m_dcaxy{-999}, m_dcaxyErr{-999}, m_dcaz{-999}, m_dcazErr{-999};
  float m_tpcSeedPt{-999}, m_tpcSeedEta{-999}, m_tpcSeedPhi{-999};  // TPC track seed (circle fit), v2
};

#endif
