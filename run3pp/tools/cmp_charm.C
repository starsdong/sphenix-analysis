// cmp_charm.C : compare two KFP 3-prong ntuples (e.g. plain 3-body vs pseudo-intermediate descriptor) by daughter track IDs
// usage: root -b -q 'cmp_charm.C+("a.root","b.root","Dplus")'
#include <TFile.h>
#include <TTree.h>

#include <algorithm>
#include <array>
#include <cstdio>
#include <map>
#include <set>
#include <string>

using Key = std::array<long, 4>;  // event, sorted track IDs

static std::map<Key, std::map<std::string, float>> load(const char *fn, int &nent, int &ndup)
{
  std::map<Key, std::map<std::string, float>> m;
  TFile f(fn);
  TTree *t = (TTree *) f.Get("DecayTree");
  nent = ndup = 0;
  if (!t) return m;
  std::string M;
  for (auto *o : *t->GetListOfBranches())
  {
    std::string n = o->GetName();
    if (n.size() > 5 && n.substr(n.size() - 5) == "_mass" && n.rfind("track_", 0) != 0)
    {
      M = n.substr(0, n.size() - 5);
      break;
    }
  }
  int ev, id[3];
  float mass, dira, chi2, pairdca[3];
  unsigned int ndof;
  t->SetBranchAddress("eventNumber", &ev);
  int nid = 0;  // daughter track-ID branches: track_N_track_ID or <intermediate>_track_N_track_ID
  for (auto *o : *t->GetListOfBranches())
  {
    std::string n = o->GetName();
    if (n.size() > 9 && n.substr(n.size() - 9) == "_track_ID" && nid < 3) t->SetBranchAddress(n.c_str(), &id[nid++]);
  }
  t->SetBranchAddress((M + "_mass").c_str(), &mass);
  t->SetBranchAddress((M + "_DIRA").c_str(), &dira);
  t->SetBranchAddress((M + "_chi2").c_str(), &chi2);
  t->SetBranchAddress((M + "_nDoF").c_str(), &ndof);
  t->SetBranchAddress("track_1_track_2_DCA", &pairdca[0]);
  t->SetBranchAddress("track_1_track_3_DCA", &pairdca[1]);
  t->SetBranchAddress("track_2_track_3_DCA", &pairdca[2]);
  nent = t->GetEntries();
  for (int e = 0; e < nent; ++e)
  {
    t->GetEntry(e);
    std::array<long, 3> s = {id[0], id[1], id[2]};
    std::sort(s.begin(), s.end());
    Key k = {ev, s[0], s[1], s[2]};
    if (m.count(k)) ++ndup;
    m[k] = {{"mass", mass}, {"dira", dira}, {"chi2ndf", ndof ? chi2 / ndof : -1}, {"maxpairdca", std::max({pairdca[0], pairdca[1], pairdca[2]})}};
  }
  printf("%s: mother %s, %d entries, %d duplicate track sets\n", fn, M.c_str(), nent, ndup);
  return m;
}

void cmp_charm(const char *fa, const char *fb, const char *tag = "")
{
  int na, nb, da, db;
  auto A = load(fa, na, da);
  auto B = load(fb, nb, db);
  int both = 0, onlyA = 0, onlyB = 0;
  double dm = 0;
  for (auto &[k, v] : A)
    if (B.count(k)) { ++both; dm = std::max(dm, (double) std::abs(v.at("mass") - B[k].at("mass"))); }
    else ++onlyA;
  for (auto &[k, v] : B)
    if (!A.count(k)) ++onlyB;
  printf("[%s] unique track sets: A %zu, B %zu, common %d (max |dm| %.4f GeV), only A %d, only B %d\n", tag, A.size(), B.size(), both, dm, onlyA, onlyB);
  int n = 0;
  for (auto &[k, v] : B)
    if (!A.count(k) && n++ < 8)
      printf("  only B: ev %ld mass %.3f DIRA %.4f chi2/ndf %.2f max pair DCA %.4f\n", k[0], v.at("mass"), v.at("dira"), v.at("chi2ndf"), v.at("maxpairdca"));
}
