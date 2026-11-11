// repro_v2.C : at-scale check that the d0-pvfix v2 production reproduces the stock candidate set.
// For every v2 file with a matching v1 (d0-pvfix v1 = stock candidate set) file: N(v2, pvfix_passOriginal==1) vs N(v1),
// plus total v2 candidates and those passing only with a refitted PV.
//   root -l -b -q 'repro_v2.C+("<v2 dir>", "<v1 dir>", "<list of v2 file names>")'
#include <TFile.h>
#include <TTree.h>

#include <cstdio>
#include <fstream>
#include <string>

void repro_v2(const char *d2, const char *d1, const char *list)
{
  std::ifstream in(list);
  std::string fn;
  long nf = 0, nmis = 0, n1 = 0, n2 = 0, n2o = 0, n2r = 0, n2k = 0, n2onlyR = 0, n2onlyK = 0, nMissing1 = 0;
  while (std::getline(in, fn))
  {
    TFile f2((std::string(d2) + "/" + fn).c_str());
    TFile f1((std::string(d1) + "/" + fn).c_str());
    auto *t2 = (TTree *) f2.Get("DecayTree");
    auto *t1 = f1.IsZombie() ? nullptr : (TTree *) f1.Get("DecayTree");
    if (!t2) continue;
    if (!t1)
    {
      ++nMissing1;
      continue;
    }
    ++nf;
    long a = t1->GetEntries(), b = t2->GetEntries(), o = t2->GetEntries("pvfix_passOriginal==1");
    n1 += a;
    n2 += b;
    n2o += o;
    n2r += t2->GetEntries("pvfix_passRefit==1");
    n2k += t2->GetEntries("pvfix_passKF==1");
    n2onlyR += t2->GetEntries("pvfix_passOriginal==0 && pvfix_passRefit==1 && pvfix_passKF==0");
    n2onlyK += t2->GetEntries("pvfix_passOriginal==0 && pvfix_passRefit==0 && pvfix_passKF==1");
    if (a != o)
    {
      if (nmis < 20) printf("MISMATCH %s: v1 %ld, v2 passOriginal %ld\n", fn.c_str(), a, o);
      ++nmis;
    }
  }
  printf("repro_v2: files compared %ld (v1 file missing for %ld) | mismatching files %ld\n", nf, nMissing1, nmis);
  printf("repro_v2: v1 candidates %ld | v2 candidates %ld (x%.3f): passOriginal %ld, passRefit %ld, passKF %ld, only refit %ld, only KF %ld\n", n1, n2,
         n1 ? (double) n2 / n1 : 0., n2o, n2r, n2k, n2onlyR, n2onlyK);
}
