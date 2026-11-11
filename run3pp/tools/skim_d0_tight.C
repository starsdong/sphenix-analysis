// skim_d0_tight.C : further D0 skim of the full-dataset skim, for local analysis
// usage: root -b -q 'skim_d0_tight.C+("in_glob_or_file", "out.root")'
#include <TChain.h>
#include <TFile.h>
#include <TTree.h>

#include <cstdio>

const char *kD0TightCut = "D0_pT>1.5";  // keeps the non-pointing DIRA 0.95-0.98 control region

void skim_d0_tight(const char *in, const char *out)
{
  TChain ch("DecayTree");
  int nf = ch.Add(in);
  Long64_t nin = ch.GetEntries();
  TFile fo(out, "RECREATE", "", 505);
  TTree *t = ch.CopyTree(kD0TightCut);
  t->SetTitle(Form("D0 tight skim: |PV DCA|>0.006, daughter pT>0.5, %s", kD0TightCut));
  Long64_t nout = t->GetEntries();
  t->Write("", TObject::kOverwrite);
  fo.Close();
  printf("skim_d0_tight: %d files, %lld -> %lld candidates (%.2f%%)\n", nf, nin, nout, nin ? 100. * nout / nin : 0.);
}
