// Diagnostic: what does SvtxTrack::get_vertex_id() point to in run3pp DST_TRKR_TRACKS?
R__LOAD_LIBRARY(libtrackbase_historic_io.so)
R__LOAD_LIBRARY(libglobalvertex_io.so)
R__LOAD_LIBRARY(libtrack_io.so)
#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/SvtxTrack.h>
#include <globalvertex/SvtxVertexMap.h>
#include <globalvertex/SvtxVertex.h>
void check_vertex_id(const char *fn, int nev = 100)
{
  TFile *f = TFile::Open(fn);
  TTree *T = (TTree *) f->Get("T");
  SvtxTrackMap *trk = nullptr;
  SvtxVertexMap *vtx = nullptr;
  T->SetBranchAddress("DST#SVTX#SvtxTrackMap", &trk);
  T->SetBranchAddress("DST#SVTX#SvtxVertexMap", &vtx);
  long n = 0, nUnset = 0, nNotInMap = 0, nInMap = 0, nSameX = 0, nInVtxList = 0, nInVtxListSameX = 0;
  long nListed = 0, nListedTrkFound = 0, nListedTrkIdBack = 0;
  TH1F *hdz = new TH1F("hdz", "", 200, -50, 50);
  TH1F *hdzSame = new TH1F("hdzSame", "", 200, -50, 50);
  for (int i = 0; i < std::min<long>(nev, T->GetEntries()); ++i)
  {
    T->GetEntry(i);
    // vertex -> track list consistency
    for (auto &[vid, v] : *vtx)
    {
      for (auto it = v->begin_tracks(); it != v->end_tracks(); ++it)
      {
        ++nListed;
        auto *t = trk->get(*it);
        if (t)
        {
          ++nListedTrkFound;
          if (t->get_vertex_id() == vid) ++nListedTrkIdBack;
        }
      }
    }
    for (auto &[id, t] : *trk)
    {
      ++n;
      unsigned int vid = t->get_vertex_id();
      if (vid == UINT_MAX) { ++nUnset; continue; }
      auto *v = vtx->get(vid);
      if (!v) { ++nNotInMap; continue; }
      ++nInMap;
      bool same = (int) v->get_beam_crossing() == (int) t->get_crossing();
      if (same) ++nSameX;
      bool listed = false;
      for (auto it = v->begin_tracks(); it != v->end_tracks(); ++it) if (*it == id) listed = true;
      if (listed) { ++nInVtxList; if (same) ++nInVtxListSameX; }
      hdz->Fill(t->get_z() - v->get_z());
      if (same && listed) hdzSame->Fill(t->get_z() - v->get_z());
    }
  }
  std::cout << "tracks " << n << " vertex_id unset " << nUnset << " not in map " << nNotInMap << " in map " << nInMap
            << " (same crossing " << nSameX << ", track listed in vertex " << nInVtxList << ", listed&same " << nInVtxListSameX << ")" << std::endl;
  std::cout << "vertex track lists: entries " << nListed << " track found " << nListedTrkFound << " track vertex_id points back " << nListedTrkIdBack << std::endl;
  std::cout << "dz(track pca - vertex) all in-map: RMS " << hdz->GetRMS() << " frac |dz|<1cm " << hdz->Integral(hdz->FindBin(-1), hdz->FindBin(1)) / hdz->GetEntries() << std::endl;
  std::cout << "dz listed&same: entries " << hdzSame->GetEntries() << " RMS " << hdzSame->GetRMS() << std::endl;
}
