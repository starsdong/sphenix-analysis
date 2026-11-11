// skim_d0.C : second-stage D0 skim of the (slimmed) prod KFP_D0 ntuples. Output has the same tree format.
// usage: root -b -q 'skim_d0.C+("files.txt","out.root")'   (files.txt: one KFP_D0_*.root per line)
// cuts: both daughters |PV_DCA| > 0.006 cm (3D, to the candidate's PV), D0_DIRA > 0.95, both daughters pT > 0.5 GeV/c
#include <TChain.h>
#include <TFile.h>
#include <TTree.h>

#include <cstdio>
#include <fstream>
#include <string>

const char *kD0SkimCut =
    "abs(track_1_PV_DCA)>0.006 && abs(track_2_PV_DCA)>0.006 && D0_DIRA>0.95 && track_1_pT>0.5 && track_2_pT>0.5";

void skim_d0(const char *list, const char *out)
{
  TChain ch("DecayTree");
  std::ifstream in(list);
  std::string f;
  int nf = 0;
  while (std::getline(in, f))
    if (!f.empty()) nf += ch.Add(f.c_str());
  Long64_t nin = ch.GetEntries();
  TFile fo(out, "RECREATE", "", 505);  // ZSTD level 5, as the slimmed production files
  TTree *t = ch.CopyTree(kD0SkimCut);
  t->SetTitle(Form("D0 skim: %s", kD0SkimCut));
  Long64_t nout = t->GetEntries();
  t->Write("", TObject::kOverwrite);
  fo.Close();
  printf("skim_d0: %d files, %lld -> %lld candidates (%.2f%%)\n", nf, nin, nout, nin ? 100. * nout / nin : 0.);
}
