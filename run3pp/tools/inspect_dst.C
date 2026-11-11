// Quick look at DST_TRKR_TRACKS content: nodes, track class, crossings, silicon clusters
R__LOAD_LIBRARY(libtrackbase_historic_io.so)
R__LOAD_LIBRARY(libglobalvertex_io.so)
R__LOAD_LIBRARY(libtrack_io.so)
#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxTrack_v5.h>
#include <globalvertex/SvtxVertexMap.h>
#include <globalvertex/SvtxVertex.h>
void inspect_dst(const char* fn, int nev = 50)
{
  TFile* f = TFile::Open(fn);
  TTree* T = (TTree*) f->Get("T");
  std::cout << "Entries " << T->GetEntries() << std::endl;
  for (auto* b : *T->GetListOfBranches()) {
    TBranch* br = (TBranch*) b;
    std::cout << "BR " << br->GetName() << " " << br->GetClassName() << " zip/ev=" << br->GetZipBytes("*")/double(T->GetEntries())/1024. << " kB" << std::endl;
  }
  SvtxTrackMap* trk = nullptr; SvtxVertexMap* vtx = nullptr;
  T->SetBranchAddress("DST#SVTX#SvtxTrackMap", &trk);
  T->SetBranchAddress("DST#SVTX#SvtxVertexMap", &vtx);
  std::map<int,int> trkX, vtxX, tpcOnlyX;
  long ntrk = 0, nsil = 0, nvtx = 0;
  for (int i = 0; i < nev; ++i) {
    T->GetEntry(i);
    if (i < 3) std::cout << "ev " << i << " ntrk " << trk->size() << " nvtx " << vtx->size() << " trk class " << (trk->size() ? trk->begin()->second->ClassName() : "") << std::endl;
    for (auto& [id, t] : *trk) {
      ++ntrk;
      auto* v5 = dynamic_cast<SvtxTrack_v5*>(t);
      int nsi = v5 ? v5->get_nmvtx_clusters() + v5->get_nintt_clusters() : -1;
      if (nsi > 0) { ++nsil; trkX[t->get_crossing()]++; }
      else tpcOnlyX[t->get_crossing()]++;
      if (i == 0 && id < 8) std::cout << "  trk " << id << " pt " << t->get_pt() << " q " << t->get_charge() << " X " << t->get_crossing() << " nmvtx " << (v5 ? v5->get_nmvtx_clusters() : -1) << " nintt " << (v5 ? v5->get_nintt_clusters() : -1) << " ntpc " << (v5 ? v5->get_ntpc_clusters() : -1) << " ntpot " << (v5 ? v5->get_ntpot_clusters() : -1) << " vtxid " << t->get_vertex_id() << " z " << t->get_z() << " chi2/ndf " << t->get_chisq() / t->get_ndf() << std::endl;
    }
    for (auto& [id, v] : *vtx) {
      ++nvtx; vtxX[v->get_beam_crossing()]++;
      if (i == 0) std::cout << "  vtx " << id << " X " << v->get_beam_crossing() << " ntrk " << v->size_tracks() << " xyz " << v->get_x() << " " << v->get_y() << " " << v->get_z() << " chi2/ndf " << v->get_chisq() << "/" << v->get_ndof() << std::endl;
    }
  }
  std::cout << "tracks/ev " << ntrk / double(nev) << " with silicon " << nsil / double(nev) << " vertices/ev " << nvtx / double(nev) << std::endl;
  auto pr = [](const char* t, std::map<int,int>& m) { std::cout << t; int k = 0; for (auto& [x, n] : m) { if (k++ < 40) std::cout << " " << x << ":" << n; } std::cout << " (" << m.size() << " distinct)" << std::endl; };
  pr("silicon-track crossings:", trkX); pr("TPC-only track crossings:", tpcOnlyX); pr("vertex crossings:", vtxX);
}
