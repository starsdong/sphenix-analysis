LD_LIBRARY_PATH now /sphenix/u/dongxs/install/lib:/cvmfs/sphenix.sdcc.bnl.gov/alma9.2-gcc-14.2.0/opt/sphenix/core/gcc/14.2.0-2f0a0/x86_64-el9/lib:/cvmfs/sphenix.sdcc.bnl.gov/alma9.2-gcc-14.2.0/opt/sphenix/core/gcc/14.2.0-2f0a0/x86_64-el9/lib64:/cvmfs/sphenix.sdcc.bnl.gov/alma9.2-gcc-14.2.0/opt/sphenix/core/binutils/2.40-acaab/x86_64-el9/lib:/cvmfs/sphenix.sdcc.bnl.gov/alma9.2-gcc-14.2.0/release/release_new/new.8/lib64:/cvmfs/sphenix.sdcc.bnl.gov/alma9.2-gcc-14.2.0/release/release_new/new.8/lib:/cvmfs/sphenix.sdcc.bnl.gov/alma9.2-gcc-14.2.0/opt/sphenix/utils/lib64:/cvmfs/sphenix.sdcc.bnl.gov/alma9.2-gcc-14.2.0/opt/sphenix/utils/lib:/cvmfs/sphenix.sdcc.bnl.gov/alma9.2-gcc-14.2.0/opt/sphenix/core/lib:/cvmfs/sphenix.sdcc.bnl.gov/alma9.2-gcc-14.2.0/opt/sphenix/core/root-6.32.06/lib:/cvmfs/sphenix.sdcc.bnl.gov/alma9.2-gcc-14.2.0/opt/sphenix/core/geant4.10.07.p04/lib64:/cvmfs/sphenix.sdcc.bnl.gov/alma9.2-gcc-14.2.0/opt/sphenix/core/lhapdf-5.9.1/lib:/usr/local/lib64:/usr/lib64
PATH now /sphenix/u/dongxs/install/bin:/cvmfs/sphenix.sdcc.bnl.gov/alma9.2-gcc-14.2.0/opt/sphenix/core/gcc/14.2.0-2f0a0/x86_64-el9/bin:/cvmfs/sphenix.sdcc.bnl.gov/alma9.2-gcc-14.2.0/opt/sphenix/core/binutils/2.40-acaab/x86_64-el9/bin:/cvmfs/sphenix.sdcc.bnl.gov/alma9.2-gcc-14.2.0/release/release_new/new.8/bin:/cvmfs/sphenix.sdcc.bnl.gov/alma9.2-gcc-14.2.0/opt/sphenix/utils/bin:/cvmfs/sphenix.sdcc.bnl.gov/alma9.2-gcc-14.2.0/opt/sphenix/core/bin:/cvmfs/sphenix.sdcc.bnl.gov/alma9.2-gcc-14.2.0/opt/sphenix/core/root-6.32.06/bin:/cvmfs/sphenix.sdcc.bnl.gov/alma9.2-gcc-14.2.0/opt/sphenix/core/geant4.10.07.p04/bin:/usr/local/bin:/usr/bin:/usr/local/sbin:/usr/sbin
// slim_kfp.C : remove redundant branches from a KFParticle "DecayTree" ntuple and recompress (ZSTD).
//
// Used by Fun4All_run3pp_KFP.C for the "prod" cut set (in-job, after se->End()), and standalone:
//   root -b -q 'slim_kfp.C+("KFP_Ks_xxx.root")'          (file is replaced in place)
//
// Kept per mother M : mass, massErr, pT, pTErr, pseudorapidity, phi, PDG_ID, all topology
//                     (x,y,z, decay length/time (+xy, errors), DIRA(+xy), FDchi2, PV_DCA*, chi2, nDoF,
//                     vertex_volume)
// Kept per daughter : pT, pTErr, pseudorapidity, phi, charge, PDG_ID, track_ID, bunch_crossing, chi2, nDoF,
//                     PV_DCA*, dEdx, detector nHits/nStates
// Dropped (recomputable from the kept ones or redundant):
//   all *_Covariance arrays (mother 21, daughters 2x21, PV 6)
//   mother   : px py pz pE p pErr theta rapidity charge
//   daughter : px py pz pE p pErr theta rapidity mass jT x y z
// Also fixes KFParticle's 6/21-byte basket size of the covariance branches (pure storage overhead).
//   event    : secondary_vertex_mass_pionPID (= mother mass for Ks), all trigger_* (64 bools; trigger info is
//              in the EventInfo tree, GL1_scaled / GL1_live)
#include <TFile.h>
#include <TObjArray.h>
#include <TRegexp.h>
#include <TString.h>
#include <TSystem.h>
#include <TTree.h>

#include <iostream>
#include <string>
#include <vector>

bool slimDecayTree(const std::string &fname)
{
  TFile *in = TFile::Open(fname.c_str());
  if (!in || in->IsZombie()) return false;
  TTree *t = (TTree *) in->Get("DecayTree");
  if (!t)
  {
    in->Close();
    return false;
  }
  // mother prefix = prefix of the first "*_mass" branch that is not a track
  std::string M;
  for (auto *o : *t->GetListOfBranches())
  {
    TString n = o->GetName();
    if (n.EndsWith("_mass") && !n.BeginsWith("track_"))
    {
      M = TString(n(0, n.Length() - 5)).Data();
      break;
    }
  }
  // intermediates (e.g. phi in Ds -> phi pi): other non-track prefixes with a "_mass" branch
  std::vector<std::string> inter;
  int ntrk = 0;
  for (auto *o : *t->GetListOfBranches())
  {
    TString n = o->GetName();
    if (n.EndsWith("_mass") && !n.BeginsWith("track_") && n != (M + "_mass").c_str()) inter.push_back(TString(n(0, n.Length() - 5)).Data());
    for (int i = 1; i <= 9; ++i)
      if (n == Form("track_%d_pT", i) && i > ntrk) ntrk = i;
  }
  std::vector<std::string> drop = {"secondary_vertex_mass_pionPID"};
  for (const char *v : {"px", "py", "pz", "pE", "p", "pErr", "theta", "rapidity", "charge"}) drop.push_back(M + "_" + v);
  for (auto &I : inter)
    for (const char *v : {"px", "py", "pz", "pE", "p", "pErr", "theta", "rapidity", "charge"}) drop.push_back(I + "_" + v);
  for (int i = 1; i <= 9; ++i)
    for (const char *v : {"px", "py", "pz", "pE", "p", "pErr", "theta", "rapidity", "mass", "jT", "x", "y", "z"})
      drop.push_back(Form("track_%d_%s", i, v));
  (void) ntrk;
  int ndrop = 0;
  for (auto *o : *t->GetListOfBranches())
  {
    TString n = o->GetName();
    bool d = n.BeginsWith("trigger_") || n.EndsWith("_Covariance");  // covariances dropped for production
    // daughters of intermediates are named <intermediate>_track_N_<var>
    for (auto &s : drop) d |= (n == s.c_str()) || (s.rfind("track_", 0) == 0 && n.EndsWith(("_" + s).c_str()));
    // d0-pvfix v2 production (KFP_SLIM_DROP_ALLREFITS=1): drop the candidate quantities w.r.t. the all-track refits;
    // their vertices (primary_vertex_refitAll_*, primary_vertex_kfAll_*) are kept
    static const bool dropAll = [] { const char *e = gSystem->Getenv("KFP_SLIM_DROP_ALLREFITS"); return e && e[0] == '1'; }();
    if (dropAll) d |= n.EndsWith("_PVrefitAll") || n.EndsWith("_PVkfAll");
    if (d)
    {
      t->SetBranchStatus(n, 0);
      ++ndrop;
    }
  }
  const int nbr = t->GetListOfBranches()->GetEntries();
  std::string tmp = fname + ".slim.tmp";
  TFile *out = TFile::Open(tmp.c_str(), "RECREATE");
  out->SetCompressionSettings(505);  // ZSTD level 5
  // KFParticle books the covariance branches with a buffer size of 6/21 *bytes* (4th argument of
  // TTree::Branch), i.e. one tiny basket per entry: ~100-350 B/candidate of pure basket overhead.
  // Clone the structure only, reset all basket sizes, then copy the entries (re-compressed, not "fast").
  TTree *c = t->CloneTree(0);
  c->SetBasketSize("*", 32000);
  c->CopyEntries(t);
  c->Write();
  const Long64_t nin = t->GetEntries(), nout = c->GetEntries();
  out->Close();
  in->Close();
  if (nin != nout)
  {
    std::cout << "slimDecayTree: entry mismatch " << nin << " vs " << nout << " for " << fname << std::endl;
    return false;
  }
  gSystem->Rename(tmp.c_str(), fname.c_str());
  std::cout << "slimDecayTree: " << fname << " entries " << nout << ", branches " << nbr << " -> " << nbr - ndrop << std::endl;
  return true;
}

void slim_kfp(const std::string &fname) { slimDecayTree(fname); }
