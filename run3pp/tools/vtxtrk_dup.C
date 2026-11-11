// vtxtrk_dup.C : are the same collisions reconstructed in consecutive DST events (streaming frames)?
// A vertex in event e is a "repeat" if event e-1 has a vertex with |dz| < 50 um, |dx|,|dy| < 50 um and the same size_tracks.
// root -l -b -q 'vtxtrk_dup.C("<audit _vtx.root>")'
#include <TFile.h>
#include <TTree.h>

#include <cmath>
#include <cstdio>
#include <map>
#include <vector>

void vtxtrk_dup(const char *audit)
{
  TFile fa(audit);
  TTree *v = (TTree *) fa.Get("vtx");
  int run, evt, cross, nAll;
  float vx, vy, vz;
  v->SetBranchAddress("run", &run);
  v->SetBranchAddress("evt", &evt);
  v->SetBranchAddress("cross", &cross);
  v->SetBranchAddress("nAll", &nAll);
  v->SetBranchAddress("vx", &vx);
  v->SetBranchAddress("vy", &vy);
  v->SetBranchAddress("vz", &vz);
  struct Vt
  {
    int cross, n;
    float x, y, z;
  };
  std::map<int, std::vector<Vt>> ev;
  for (long i = 0; i < v->GetEntries(); ++i)
  {
    v->GetEntry(i);
    ev[evt].push_back({cross, nAll, vx, vy, vz});
  }
  long nv = 0, rep = 0, repSame = 0, nvWithPrev = 0;
  std::map<int, int> shift;
  for (auto &[e, vs] : ev)
  {
    auto p = ev.find(e - 1);
    for (auto &a : vs)
    {
      ++nv;
      if (p == ev.end()) continue;
      ++nvWithPrev;
      for (auto &b : p->second)
        if (std::fabs(a.z - b.z) < 0.005 && std::fabs(a.x - b.x) < 0.005 && std::fabs(a.y - b.y) < 0.005)
        {
          ++rep;
          repSame += (a.n == b.n);
          shift[a.cross - b.cross]++;
          break;
        }
    }
  }
  // same check within one event, different crossings
  long within = 0;
  for (auto &[e, vs] : ev)
    for (size_t i = 0; i < vs.size(); ++i)
      for (size_t j = 0; j < i; ++j)
        if (vs[i].cross != vs[j].cross && std::fabs(vs[i].z - vs[j].z) < 0.005 && std::fabs(vs[i].x - vs[j].x) < 0.005) ++within;
  printf("vertices %ld (events %zu), with a previous event %ld; repeated from the previous event %ld (%.1f%%), of those same size_tracks %ld\n", nv,
         ev.size(), nvWithPrev, rep, 100. * rep / nvWithPrev, repSame);
  printf("pairs of vertices within one event, different crossings, |dx|,|dz| < 50 um: %ld\n", within);
  printf("crossing shift (this - previous) of repeated vertices:");
  int k = 0;
  for (auto &[s, n] : shift)
    if (n > 2 && k++ < 15) printf(" %d:%d", s, n);
  printf("\n");
}
