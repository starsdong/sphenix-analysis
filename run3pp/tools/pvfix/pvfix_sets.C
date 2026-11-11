// pvfix_sets.C : compare candidate sets of the smoke tests (key = event, track_1 ID, track_2 ID)
#include <TFile.h>
#include <TTree.h>

#include <array>
#include <cstdio>
#include <map>
#include <set>
#include <tuple>

namespace
{
  using K = std::tuple<int, int, int>;
  using P = std::tuple<int, int>;  // event, unordered pair hash
  std::map<K, std::array<int, 3>> load(const char *f)
  {
    TFile F(f);
    TTree *t = (TTree *) F.Get("DecayTree");
    int e, a, b, po = 1, pr = 0, pk = 0;
    t->SetBranchAddress("eventNumber", &e);
    t->SetBranchAddress("track_1_track_ID", &a);
    t->SetBranchAddress("track_2_track_ID", &b);
    bool has = t->GetBranch("pvfix_passOriginal");
    if (has)
    {
      t->SetBranchAddress("pvfix_passOriginal", &po);
      t->SetBranchAddress("pvfix_passRefit", &pr);
      t->SetBranchAddress("pvfix_passKF", &pk);
    }
    std::map<K, std::array<int, 3>> m;
    for (long i = 0; i < t->GetEntries(); ++i)
    {
      t->GetEntry(i);
      m[{e, a, b}] = {po, pr, pk};
    }
    return m;
  }
  std::set<std::tuple<int, int, int>> pairs(const std::map<K, std::array<int, 3>> &m)
  {
    std::set<std::tuple<int, int, int>> s;
    for (auto &[k, v] : m) s.insert({std::get<0>(k), std::min(std::get<1>(k), std::get<2>(k)), std::max(std::get<1>(k), std::get<2>(k))});
    return s;
  }
}  // namespace

void pvfix_sets(const char *stock, const char *A, const char *B)
{
  auto s = load(stock), a = load(A), b = load(B);
  for (auto [name, m] : {std::pair{"A", &a}, std::pair{"B", &b}})
  {
    long po = 0, pr = 0, pk = 0, onlyR = 0, onlyK = 0, poInStock = 0, stockMissing = 0, stockNotPO = 0;
    for (auto &[k, v] : *m)
    {
      po += v[0];
      pr += v[1];
      pk += v[2];
      if (!v[0] && v[1] && !v[2]) ++onlyR;
      if (!v[0] && !v[1] && v[2]) ++onlyK;
      if (v[0] && s.count(k)) ++poInStock;
    }
    for (auto &[k, v] : s)
    {
      auto it = m->find(k);
      if (it == m->end()) ++stockMissing;
      else if (!it->second[0]) ++stockNotPO;
    }
    printf("%s: %zu cand | passOriginal %ld (in stock %ld) | passRefit %ld | passKF %ld | only refit %ld, only KF %ld | stock cand missing %ld, present but passOriginal=0 %ld\n",
           name, m->size(), po, poInStock, pr, pk, onlyR, onlyK, stockMissing, stockNotPO);
  }
  // A vs B
  long aNotB = 0, bNotA = 0, pairANotB = 0, swapped = 0;
  auto pa = pairs(a), pb = pairs(b);
  for (auto &[k, v] : a)
    if (!b.count(k)) ++aNotB;
  for (auto &[k, v] : b)
    if (!a.count(k)) ++bNotA;
  for (auto &p : pa)
    if (!pb.count(p)) ++pairANotB;
  for (auto &[k, v] : a)
    if (!b.count(k) && b.count({std::get<0>(k), std::get<2>(k), std::get<1>(k)})) ++swapped;
  printf("A not in B %ld (as track pair: %ld; same pair with K/pi swapped in B: %ld) | B not in A %ld\n", aNotB, pairANotB, swapped, bNotA);
}
