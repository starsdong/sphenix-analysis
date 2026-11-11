// branch_sizes.C : compressed size per branch group of a KFP DecayTree, in bytes per candidate.
// usage: root -b -q 'branch_sizes.C+("KFP_Ks_xxx.root")'
#include <TBranch.h>
#include <TFile.h>
#include <TString.h>
#include <TTree.h>

#include <cstdio>
#include <map>
#include <string>

void branch_sizes(const char *fn)
{
  TFile f(fn);
  TTree *t = (TTree *) f.Get("DecayTree");
  double n = t->GetEntries(), tot = 0;
  std::map<std::string, double> g;
  std::map<std::string, int> nb;
  for (auto *o : *t->GetListOfBranches())
  {
    TBranch *b = (TBranch *) o;
    TString s = b->GetName();
    std::string k;
    if (s.Contains("Covariance"))
      k = s.BeginsWith("track_") ? "daughter covariance (2x21)" : (s.BeginsWith("primary") ? "PV covariance (6)" : "mother covariance (21)");
    else if (s.BeginsWith("trigger_"))
      k = "trigger bits";
    else if (s.BeginsWith("track_1_track_2"))
      k = "pair DCA";
    else if (s.BeginsWith("track_") && (s.Contains("nHits") || s.Contains("nStates") || s.Contains("dEdx")))
      k = "daughter detector/dEdx";
    else if (s.BeginsWith("track_"))
      k = "daughter kinematics+quality+DCA";
    else if (s.BeginsWith("primary_vertex"))
      k = "PV";
    else if (s.Contains("runNumber") || s.Contains("eventNumber") || s.Contains("BCO") || s.Contains("OfBC") || s.Contains("OfVertex") || s.Contains("Multiplicity") || s.Contains("secondary"))
      k = "event";
    else
      k = "mother kinematics+topology";
    g[k] += b->GetZipBytes("*");
    nb[k]++;
    tot += b->GetZipBytes("*");
  }
  printf("%s: %.0f candidates, %d branches, %.0f B/cand (branch zip bytes)\n", fn, n, t->GetListOfBranches()->GetEntries(), tot / n);
  for (auto &[k, v] : g) printf("  %-34s %3d branches %7.1f B/cand  %5.1f%%\n", k.c_str(), nb[k], v / n, 100 * v / tot);
}
