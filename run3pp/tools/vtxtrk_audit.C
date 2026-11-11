// vtxtrk_audit.C : what are the tracks attached to the SvtxVertex objects (ana581 DST_TRKR_TRACKS)?
// run on SDCC with the ana.581 environment:
//   root -l -b -q 'vtxtrk_audit.C+("<dst list file>", "<out prefix>", <max events per file>)'
// For every vertex: size_tracks, and counts after cumulative cleaning steps; per crossing: all tracks in the
// crossing; per track (first 2000 vertices): details. Writes <prefix>.txt (summary) and <prefix>_vtx.root (tree).
#include <ffaobjects/EventHeader.h>
#include <globalvertex/SvtxVertex.h>
#include <globalvertex/SvtxVertexMap.h>
#include <trackbase/TrkrDefs.h>
#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/TrackSeed.h>
#include <trackbase_historic/TrackSeedContainer.h>

#include <TFile.h>
#include <TH1.h>
#include <TSystem.h>
#include <TTree.h>

#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <vector>

R__LOAD_LIBRARY(libtrackbase_historic_io.so)
R__LOAD_LIBRARY(libglobalvertex_io.so)
R__LOAD_LIBRARY(libtrack_io.so)
R__LOAD_LIBRARY(libffaobjects.so)

namespace
{
  int nclus(const TrackSeed *s, TrkrDefs::TrkrId id)
  {
    if (!s) return 0;
    int n = 0;
    for (auto it = s->begin_cluster_keys(); it != s->end_cluster_keys(); ++it)
      if (TrkrDefs::getTrkrId(*it) == id) ++n;
    return n;
  }
  // 3D distance of closest approach of the straight-line extrapolation of the track state to the vertex
  double dca3(const SvtxTrack *t, const SvtxVertex *v)
  {
    double dx = t->get_x() - v->get_x(), dy = t->get_y() - v->get_y(), dz = t->get_z() - v->get_z();
    double px = t->get_px(), py = t->get_py(), pz = t->get_pz(), p2 = px * px + py * py + pz * pz;
    double s = (dx * px + dy * py + dz * pz) / p2;
    double ex = dx - s * px, ey = dy - s * py, ez = dz - s * pz;
    return std::sqrt(ex * ex + ey * ey + ez * ez);
  }
}  // namespace

void vtxtrk_audit(const char *list, const char *prefix, int maxEv = 200)
{
  std::ifstream in(list);
  std::string fn;
  TFile out(Form("%s_vtx.root", prefix), "RECREATE");
  TTree tv("vtx", "per vertex");
  int run, evt, vid, cross, nAll, nPt, nEta, nQual, nSil, nTpc, nInttHit, nDcaOk, nUniqSi, nUniqTpc, nClean, nCross, nVtxInCross, nCrossClean;
  float vx, vy, vz;
  tv.Branch("run", &run);
  tv.Branch("evt", &evt);
  tv.Branch("vid", &vid);
  tv.Branch("cross", &cross);
  tv.Branch("vx", &vx);
  tv.Branch("vy", &vy);
  tv.Branch("vz", &vz);
  tv.Branch("nAll", &nAll);            // size_tracks()
  tv.Branch("nPt", &nPt);              // + pT > 0.2
  tv.Branch("nEta", &nEta);            // + |eta| < 1.1
  tv.Branch("nQual", &nQual);          // + chi2/ndf < 100 and quality finite
  tv.Branch("nSil", &nSil);            // + >= 3 MVTX and >= 1 INTT clusters
  tv.Branch("nTpc", &nTpc);            // + >= 20 TPC clusters
  tv.Branch("nDcaOk", &nDcaOk);        // + 3D DCA to vertex < 0.1 cm
  tv.Branch("nUniqSi", &nUniqSi);      // + not sharing its silicon seed with an earlier track in the vertex
  tv.Branch("nUniqTpc", &nUniqTpc);    // + not sharing its TPC seed with an earlier track
  tv.Branch("nClean", &nClean);        // = nUniqTpc
  tv.Branch("nInttHit", &nInttHit);    // tracks in the vertex with >= 1 INTT cluster (crossing from INTT timing)
  tv.Branch("nCross", &nCross);        // all tracks of the SvtxTrackMap with this crossing
  tv.Branch("nCrossClean", &nCrossClean);
  tv.Branch("nVtxInCross", &nVtxInCross);
  TTree tt("trk", "tracks of the first vertices");
  float pt, eta, qual, dca;
  int nm, ni, nt, sid, tid, tcross, inVtx, tvid;
  tt.Branch("run", &run);
  tt.Branch("evt", &evt);
  tt.Branch("vid", &tvid);
  tt.Branch("pt", &pt);
  tt.Branch("eta", &eta);
  tt.Branch("qual", &qual);
  tt.Branch("dca", &dca);
  tt.Branch("nm", &nm);
  tt.Branch("ni", &ni);
  tt.Branch("nt", &nt);
  tt.Branch("sid", &sid);
  tt.Branch("tid", &tid);
  tt.Branch("cross", &tcross);
  long nev = 0, nvtx = 0, ntrkSeedMissing = 0;
  while (std::getline(in, fn))
  {
    if (fn.empty()) continue;
    TFile *f = TFile::Open(fn.c_str());
    if (!f || f->IsZombie()) continue;
    TTree *T = (TTree *) f->Get("T");
    SvtxTrackMap *tm = nullptr;
    SvtxVertexMap *vm = nullptr;
    TrackSeedContainer *sic = nullptr, *tpc = nullptr;
    EventHeader *eh = nullptr;
    T->SetBranchAddress("DST#SVTX#SvtxTrackMap", &tm);
    T->SetBranchAddress("DST#SVTX#SvtxVertexMap", &vm);
    T->SetBranchAddress("DST#SVTX#SiliconTrackSeedContainer", &sic);
    T->SetBranchAddress("DST#SVTX#TpcTrackSeedContainer", &tpc);
    T->SetBranchAddress("DST#EventHeader", &eh);
    long ne = std::min<long>(T->GetEntries(), maxEv);
    for (long ie = 0; ie < ne; ++ie)
    {
      T->GetEntry(ie);
      if (!tm || !vm) continue;
      ++nev;
      run = eh ? eh->get_RunNumber() : -1;
      evt = eh ? eh->get_EvtSequence() : (int) ie;
      // per crossing: all tracks and clean tracks
      auto clean = [&](const SvtxTrack *t, const SvtxVertex *v, int &m, int &i, int &tp) {
        const TrackSeed *ss = (sic && t->get_silicon_seed_id() != USHRT_MAX) ? sic->get(t->get_silicon_seed_id()) : t->get_silicon_seed();
        const TrackSeed *ts = (tpc && t->get_tpc_seed_id() != USHRT_MAX) ? tpc->get(t->get_tpc_seed_id()) : t->get_tpc_seed();
        if (!ss && !ts) ++ntrkSeedMissing;
        m = nclus(ss, TrkrDefs::mvtxId);
        i = nclus(ss, TrkrDefs::inttId);
        tp = nclus(ts, TrkrDefs::tpcId);
        int lvl = 0;
        if (t->get_pt() > 0.2) lvl = 1; else return lvl;
        if (std::fabs(t->get_eta()) < 1.1) lvl = 2; else return lvl;
        if (std::isfinite(t->get_quality()) && t->get_quality() < 100) lvl = 3; else return lvl;
        if (m >= 3 && i >= 1) lvl = 4; else return lvl;
        if (tp >= 20) lvl = 5; else return lvl;
        if (!v || dca3(t, v) < 0.1) lvl = 6;
        return lvl;
      };
      std::map<int, int> crossAll, crossClean, vtxPerCross;
      for (auto &[k, t] : *tm)
      {
        int m, i, tp;
        crossAll[t->get_crossing()]++;
        if (clean(t, nullptr, m, i, tp) >= 5) crossClean[t->get_crossing()]++;
      }
      for (auto &[k, v] : *vm) vtxPerCross[v->get_beam_crossing()]++;
      for (auto &[k, v] : *vm)
      {
        ++nvtx;
        vid = k;
        cross = v->get_beam_crossing();
        vx = v->get_x();
        vy = v->get_y();
        vz = v->get_z();
        nAll = v->size_tracks();
        nPt = nEta = nQual = nSil = nTpc = nDcaOk = nUniqSi = nUniqTpc = nInttHit = 0;
        std::set<int> sids, tids;
        for (auto it = v->begin_tracks(); it != v->end_tracks(); ++it)
        {
          SvtxTrack *t = tm->get(*it);
          if (!t) continue;
          int m, i, tp;
          int lvl = clean(t, v, m, i, tp);
          if (i >= 1) ++nInttHit;
          nPt += lvl >= 1;
          nEta += lvl >= 2;
          nQual += lvl >= 3;
          nSil += lvl >= 4;
          nTpc += lvl >= 5;
          nDcaOk += lvl >= 6;
          if (lvl >= 6)
          {
            int s = t->get_silicon_seed_id(), q = t->get_tpc_seed_id();
            if (!sids.count(s))
            {
              ++nUniqSi;
              if (!tids.count(q)) ++nUniqTpc;
            }
            sids.insert(s);
            tids.insert(q);
          }
          if (nvtx <= 2000)
          {
            tvid = k;
            pt = t->get_pt();
            eta = t->get_eta();
            qual = t->get_quality();
            dca = dca3(t, v);
            nm = m;
            ni = i;
            nt = tp;
            sid = t->get_silicon_seed_id();
            tid = t->get_tpc_seed_id();
            tcross = t->get_crossing();
            tt.Fill();
          }
        }
        nClean = nUniqTpc;
        nCross = crossAll[cross];
        nCrossClean = crossClean[cross];
        nVtxInCross = vtxPerCross[cross];
        tv.Fill();
      }
    }
    f->Close();
  }
  out.Write();
  std::ofstream s(Form("%s.txt", prefix));
  s << "events " << nev << " vertices " << nvtx << " tracks without any seed " << ntrkSeedMissing << "\n";
  for (const char *b : {"nAll", "nPt", "nEta", "nQual", "nSil", "nTpc", "nDcaOk", "nUniqSi", "nUniqTpc", "nInttHit", "nCross", "nCrossClean", "nVtxInCross"})
  {
    tv.Draw(Form("%s>>h(400,0,400)", b), "", "goff");
    auto *h = (TH1 *) gDirectory->Get("h");
    double q[3], p[3] = {0.5, 0.9, 0.99};
    h->GetQuantiles(3, q, p);
    s << Form("%-12s mean %6.2f  median %4.0f  90%% %4.0f  99%% %4.0f  frac>=20 %5.2f%%\n", b, h->GetMean(), q[0], q[1], q[2],
              100. * h->Integral(21, 401) / h->Integral(0, 401));
    delete h;
  }
  s.close();
  std::cout << "done " << nev << " events " << nvtx << " vertices" << std::endl;
}
