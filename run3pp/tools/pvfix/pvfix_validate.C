// pvfix_validate.C : validate the d0-pvfix library against the stock one on the same events.
//   root -l -b -q 'pvfix_validate.C+("stock.root","closure.root","pvfix.root")'
// Candidates are matched by (run, event, track_1_track_ID, track_2_track_ID).
// Checks: (1) same candidates in all three files; (2) every branch common to stock and pvfix is identical except
// primary_vertex_ID; (3) fix 1: pvfix primary_vertex_crossing == stock primary_vertex_ID; (4) closure file:
// *_PVrefit quantities == the original ones.
#include <TFile.h>
#include <TLeaf.h>
#include <TTree.h>

#include <cmath>
#include <cstdio>
#include <map>
#include <set>
#include <string>
#include <tuple>
#include <vector>

namespace
{
  using Key = std::tuple<int, int, int, int>;
  std::map<Key, long> index(TTree *t)
  {
    int r, e, a, b;
    t->SetBranchStatus("*", 0);
    for (auto n : {"runNumber", "eventNumber", "track_1_track_ID", "track_2_track_ID"}) t->SetBranchStatus(n, 1);
    t->SetBranchAddress("runNumber", &r);
    t->SetBranchAddress("eventNumber", &e);
    t->SetBranchAddress("track_1_track_ID", &a);
    t->SetBranchAddress("track_2_track_ID", &b);
    std::map<Key, long> m;
    for (long i = 0; i < t->GetEntries(); ++i)
    {
      t->GetEntry(i);
      m[{r, e, a, b}] = i;
    }
    t->ResetBranchAddresses();
    t->SetBranchStatus("*", 1);
    return m;
  }
  double val(TTree *t, const char *br, long i, int k = 0)
  {
    t->GetEntry(i);
    TLeaf *l = t->GetLeaf(br);
    return l->GetValue(k);
  }
}  // namespace

void pvfix_validate(const char *fs, const char *fc, const char *fp)
{
  TFile Fs(fs), Fc(fc), Fp(fp);
  TTree *s = (TTree *) Fs.Get("DecayTree"), *c = (TTree *) Fc.Get("DecayTree"), *p = (TTree *) Fp.Get("DecayTree");
  auto ms = index(s), mc = index(c), mp = index(p);
  long common = 0;
  for (auto &[k, i] : ms) common += mp.count(k) && mc.count(k);
  printf("candidates: stock %zu closure %zu pvfix %zu, in all three %ld\n", ms.size(), mc.size(), mp.size(), common);
  // (2) common branches
  std::vector<std::string> br;
  std::set<std::string> pb;
  for (auto *o : *p->GetListOfBranches()) pb.insert(o->GetName());
  for (auto *o : *s->GetListOfBranches())
    if (pb.count(o->GetName())) br.push_back(o->GetName());
  printf("branches: stock %d pvfix %d common %zu; new in pvfix:", s->GetListOfBranches()->GetEntries(), p->GetListOfBranches()->GetEntries(), br.size());
  std::set<std::string> sb;
  for (auto *o : *s->GetListOfBranches()) sb.insert(o->GetName());
  for (auto &n : pb)
    if (!sb.count(n)) printf(" %s", n.c_str());
  printf("\n");
  int nDiffBr = 0;
  for (auto &b : br)
  {
    int len = s->GetLeaf(b.c_str()) ? s->GetLeaf(b.c_str())->GetLen() : 1;
    long nd = 0, n = 0;
    for (auto &[k, i] : ms)
    {
      auto it = mp.find(k);
      if (it == mp.end()) continue;
      for (int j = 0; j < len; ++j)
      {
        double a = val(s, b.c_str(), i, j), q = val(p, b.c_str(), it->second, j);
        ++n;
        if (!(a == q || (std::isnan(a) && std::isnan(q)))) ++nd;
      }
    }
    if (nd)
    {
      ++nDiffBr;
      printf("  differs: %-30s %ld of %ld values\n", b.c_str(), nd, n);
    }
  }
  printf("common branches with any difference: %d\n", nDiffBr);
  // (3) fix 1
  long okCross = 0, n3 = 0, idSame = 0;
  for (auto &[k, i] : ms)
  {
    auto it = mp.find(k);
    if (it == mp.end()) continue;
    ++n3;
    okCross += val(s, "primary_vertex_ID", i) == val(p, "primary_vertex_crossing", it->second);
    idSame += val(s, "primary_vertex_ID", i) == val(p, "primary_vertex_ID", it->second);
  }
  printf("fix 1: pvfix primary_vertex_crossing == stock primary_vertex_ID for %ld of %ld; new ID == old ID for %ld\n", okCross, n3, idSame);
  // (4) closure
  std::string M = "D0";
  std::vector<std::pair<std::string, std::string>> pairs = {
      {M + "_DIRA", M + "_DIRA_PVrefit"},       {M + "_DIRA_xy", M + "_DIRA_xy_PVrefit"},       {M + "_FDchi2", M + "_FDchi2_PVrefit"},
      {M + "_PV_DCA", M + "_PV_DCA_PVrefit"},   {M + "_PV_DCA_xy", M + "_PV_DCA_xy_PVrefit"},   {M + "_decayLength", M + "_decayLength_PVrefit"},
      {M + "_decayLengthErr", M + "_decayLengthErr_PVrefit"}, {M + "_decayLength_xy", M + "_decayLength_xy_PVrefit"},
      {"track_1_PV_DCA", "track_1_PV_DCA_PVrefit"}, {"track_2_PV_DCA", "track_2_PV_DCA_PVrefit"}, {"track_1_PV_DCA_xy", "track_1_PV_DCA_xy_PVrefit"},
      {"track_2_PV_DCA_xy", "track_2_PV_DCA_xy_PVrefit"}, {"track_1_PV_DCA_sig", "track_1_PV_DCA_sig_PVrefit"},
      {"primary_vertex_x", "primary_vertex_refit_x"}, {"primary_vertex_z", "primary_vertex_refit_z"}};
  printf("closure (KFP_PVREFIT_CLOSURE=1): original vs *_PVrefit\n");
  for (auto &[o, r] : pairs)
  {
    if (!c->GetLeaf(o.c_str()) || !c->GetLeaf(r.c_str()))
    {
      printf("  %-28s missing\n", o.c_str());
      continue;
    }
    long n = 0, same = 0;
    double maxrel = 0;
    for (auto &[k, i] : mc)
    {
      double a = val(c, o.c_str(), i), q = val(c, r.c_str(), i);
      ++n;
      if (a == q) ++same;
      else maxrel = std::max(maxrel, std::fabs(a - q) / std::max(1e-12, std::fabs(a)));
    }
    printf("  %-28s identical %ld of %ld, max rel. diff %.2e\n", o.c_str(), same, n, maxrel);
  }
}
