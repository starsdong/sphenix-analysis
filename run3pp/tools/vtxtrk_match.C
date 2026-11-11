// vtxtrk_match.C : match KFParticle D0 candidates to the audited vertices by (run, eventNumber, primary_vertex_ID)
// root -l -b -q 'vtxtrk_match.C("<audit _vtx.root>", "<KFP D0 file>")'
#include <TFile.h>
#include <TTree.h>

#include <cstdio>
#include <map>
#include <tuple>

void vtxtrk_match(const char *audit, const char *kfp)
{
  TFile fa(audit);
  TTree *v = (TTree *) fa.Get("vtx");
  int run, evt, vid, nAll, nClean;
  v->SetBranchAddress("run", &run);
  v->SetBranchAddress("evt", &evt);
  v->SetBranchAddress("vid", &vid);
  v->SetBranchAddress("nAll", &nAll);
  v->SetBranchAddress("nClean", &nClean);
  std::map<std::tuple<int, int, int>, std::pair<int, int>> m;
  int emin = 1 << 30, emax = -1;
  for (long i = 0; i < v->GetEntries(); ++i)
  {
    v->GetEntry(i);
    m[{run, evt, vid}] = {nAll, nClean};
    emin = std::min(emin, evt);
    emax = std::max(emax, evt);
  }
  TFile fk(kfp);
  TTree *t = (TTree *) fk.Get("DecayTree");
  int r, e, p, nv;
  t->SetBranchAddress("runNumber", &r);
  t->SetBranchAddress("eventNumber", &e);
  t->SetBranchAddress("primary_vertex_ID", &p);
  t->SetBranchAddress("nTracksOfVertex", &nv);
  long inRange = 0, matched = 0, same = 0;
  double sA = 0, sC = 0, sNv = 0;
  std::map<std::tuple<int, int, int>, int> candPerVtx;
  for (long i = 0; i < t->GetEntries(); ++i)
  {
    t->GetEntry(i);
    if (e < emin || e > emax) continue;
    ++inRange;
    auto it = m.find({r, e, p});
    if (it == m.end()) continue;
    ++matched;
    same += (it->second.first == nv);
    sA += it->second.first;
    sC += it->second.second;
    sNv += nv;
    candPerVtx[{r, e, p}]++;
  }
  // vertex-weighted means
  double vA = 0, vC = 0, vAw = 0, w = 0;
  for (auto &[k, x] : m)
  {
    vA += x.first;
    vC += x.second;
    double ww = x.first * (x.first - 1.0);  // ~ number of track pairs
    vAw += ww * x.first;
    w += ww;
  }
  printf("audit: %zu vertices in events %d-%d\n", m.size(), emin, emax);
  printf("KFP candidates in that event range %ld, matched to an audited vertex %ld (%.1f%%), nTracksOfVertex == size_tracks for %.1f%%\n", inRange, matched,
         100. * matched / std::max(1L, inRange), 100. * same / std::max(1L, matched));
  printf("per candidate: <nTracksOfVertex> %.2f, <size_tracks> %.2f, <clean> %.2f ; candidates per matched vertex %.2f over %zu vertices\n", sNv / matched,
         sA / matched, sC / matched, (double) matched / candPerVtx.size(), candPerVtx.size());
  printf("per vertex (all audited): <size_tracks> %.2f, <clean> %.2f ; pair-weighted <size_tracks> %.2f\n", vA / m.size(), vC / m.size(), vAw / w);
}
