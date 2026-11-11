// trkstate_check.C : are the DST track states at the PCA to their vertex (PHActsVertexPropagator)?
// root -l -b -q 'trkstate_check.C+("<DST>", 200)'
#include <globalvertex/SvtxVertex.h>
#include <globalvertex/SvtxVertexMap.h>
#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxTrackMap.h>

#include <TFile.h>
#include <TH1D.h>
#include <TTree.h>

#include <cmath>
#include <cstdio>

R__LOAD_LIBRARY(libtrackbase_historic_io.so)
R__LOAD_LIBRARY(libglobalvertex_io.so)

void trkstate_check(const char *dst, int maxEv = 200)
{
  TFile *f = TFile::Open(dst);
  TTree *T = (TTree *) f->Get("T");
  SvtxTrackMap *tm = nullptr;
  SvtxVertexMap *vm = nullptr;
  T->SetBranchAddress("DST#SVTX#SvtxTrackMap", &tm);
  T->SetBranchAddress("DST#SVTX#SvtxVertexMap", &vm);
  TH1D hd("hd", "", 2000, 0, 2), hl("hl", "", 2000, -1, 1), hr("hr", "", 2000, 0, 20), hdo("hdo", "", 2000, 0, 2);
  long nT = 0, nVid = 0, nVidMatch = 0, nIn = 0;
  for (long ie = 0; ie < std::min<long>(T->GetEntries(), maxEv); ++ie)
  {
    T->GetEntry(ie);
    for (auto &[vk, v] : *vm)
      for (auto it = v->begin_tracks(); it != v->end_tracks(); ++it)
      {
        SvtxTrack *t = tm->get(*it);
        if (!t) continue;
        ++nIn;
        nVidMatch += (t->get_vertex_id() == vk);
        double r[3] = {t->get_x() - v->get_x(), t->get_y() - v->get_y(), t->get_z() - v->get_z()};
        double p = t->get_p(), u[3] = {t->get_px() / p, t->get_py() / p, t->get_pz() / p};
        double along = r[0] * u[0] + r[1] * u[1] + r[2] * u[2];
        double dist = std::sqrt(r[0] * r[0] + r[1] * r[1] + r[2] * r[2]);
        hd.Fill(dist);                       // |state point - vertex| (cm)
        hl.Fill(along);                      // component along the momentum: 0 if the state is at the PCA to the vertex
        hr.Fill(std::hypot(t->get_x(), t->get_y()));  // transverse radius of the state point
      }
    for (auto &[tk, t] : *tm)
    {
      ++nT;
      if (t->get_vertex_id() != UINT_MAX) ++nVid;
    }
  }
  auto q = [](TH1D &h, double p) { double x; h.GetQuantiles(1, &x, &p); return x; };
  printf("tracks %ld, with a vertex id %ld; vertex tracks %ld, whose vertex_id == their vertex %ld (%.1f%%)\n", nT, nVid, nIn, nVidMatch, 100. * nVidMatch / nIn);
  printf("vertex tracks: |state - vertex| median %.1f um, 90%% %.1f um; component along p: median |.| %.2f um, 99%% within %.2f um; state radius median %.3f cm\n",
         1e4 * q(hd, 0.5), 1e4 * q(hd, 0.9), 1e4 * std::fabs(q(hl, 0.5)), 1e4 * std::max(std::fabs(q(hl, 0.005)), std::fabs(q(hl, 0.995))), q(hr, 0.5));
}
