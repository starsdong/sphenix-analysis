// fill_v0_mass.C -- reduce KFParticle Ks / Lambda ntuples (loosecov production) to invariant-mass vs pT histograms
// for the DDSim momentum calibration. Only histograms are kept (the ntuples are deleted by the job wrapper).
//   root -l -b -q 'fill_v0_mass.C("<kfp dir>","<out.root>")'
// Cuts (no covariance or significance requirement):
//   all histograms : |PV z| < 10 cm, |eta_daughter| < 1.1
//   "qa"           : + DIRA > 0.995, L_xy > 0.5 cm, DCA12 < 0.2 cm, |DCA_xy(daughter)| > 0.02 cm (run3pp V0 QA cuts)
// Lambda: particle / antiparticle from Lambda0_PDG_ID; the proton is the daughter with mass > 0.5.
#include <TFile.h>
#include <TTree.h>
#include <TH2F.h>
#include <TH1D.h>
#include <TParameter.h>
#include <TSystem.h>
#include <TSystemDirectory.h>
#include <TList.h>
#include <TChain.h>
#include <cmath>
#include <iostream>
#include <string>

namespace fv {
long fill(TChain& t, const std::string& M, TH2F* hQa, TH2F* hLoose, TH2F* hQaAnti, TH2F* hLooseAnti,
          TH1D* hDca12, TH1D* hRxy)
{
  float mass, x, y, dira, pt, lxy, d12, dxy1, dxy2, pvz, eta1, eta2;
  int pdg = 0;
  t.SetBranchStatus("*", 0);
  auto B = [&](const std::string& n, void* a) { t.SetBranchStatus(n.c_str(), 1); t.SetBranchAddress(n.c_str(), a); };
  B(M + "_mass", &mass); B(M + "_x", &x); B(M + "_y", &y); B(M + "_DIRA", &dira); B(M + "_pT", &pt);
  B(M + "_decayLength_xy", &lxy); B(M + "_PDG_ID", &pdg);
  B("track_1_track_2_DCA", &d12); B("track_1_PV_DCA_xy", &dxy1); B("track_2_PV_DCA_xy", &dxy2);
  B("primary_vertex_z", &pvz); B("track_1_pseudorapidity", &eta1); B("track_2_pseudorapidity", &eta2);
  long n = 0;
  for (Long64_t i = 0; i < t.GetEntries(); ++i) {
    t.GetEntry(i);
    if (std::fabs(pvz) > 10 || std::fabs(eta1) > 1.1 || std::fabs(eta2) > 1.1) continue;
    const bool anti = hQaAnti && pdg < 0;
    (anti ? hLooseAnti : hLoose)->Fill(pt, mass);
    if (dira < 0.995 || lxy < 0.5 || d12 > 0.2 || std::fabs(dxy1) < 0.02 || std::fabs(dxy2) < 0.02) continue;
    (anti ? hQaAnti : hQa)->Fill(pt, mass);
    hDca12->Fill(d12); hRxy->Fill(std::hypot(x, y));
    ++n;
  }
  return n;
}
}  // namespace fv

void fill_v0_mass(const char* kfpDir, const char* outFile)
{
  TChain tk("DecayTree"), tl("DecayTree"), te("EventInfo");
  tk.Add(Form("%s/KFP_Ks_*.root", kfpDir));
  tl.Add(Form("%s/KFP_Lambda_*.root", kfpDir));
  te.Add(Form("%s/KFP_EventInfo_*.root", kfpDir));
  TFile out(outFile, "RECREATE");
  TH2F* ksQa = new TH2F("hKs_qa", "K0S, QA cuts;p_{T} (GeV/c);M (GeV/c^{2})", 100, 0, 5, 200, 0.40, 0.60);
  TH2F* ksLo = new TH2F("hKs_loose", "K0S, |PVz|<10 |eta|<1.1 only;p_{T} (GeV/c);M (GeV/c^{2})", 100, 0, 5, 200, 0.40, 0.60);
  TH2F* laQa = new TH2F("hLa_qa", "Lambda, QA cuts;p_{T};M", 100, 0, 5, 320, 1.08, 1.16);
  TH2F* laLo = new TH2F("hLa_loose", "Lambda, loose;p_{T};M", 100, 0, 5, 320, 1.08, 1.16);
  TH2F* lbQa = new TH2F("hLaBar_qa", "anti-Lambda, QA cuts;p_{T};M", 100, 0, 5, 320, 1.08, 1.16);
  TH2F* lbLo = new TH2F("hLaBar_loose", "anti-Lambda, loose;p_{T};M", 100, 0, 5, 320, 1.08, 1.16);
  TH1D* d12k = new TH1D("hKs_dca12_qa", "K0S QA, all masses;DCA12 (cm)", 200, 0, 0.2);
  TH1D* rk = new TH1D("hKs_rxy_qa", "K0S QA, all masses;R_{xy} (cm)", 200, 0, 20);
  TH1D* d12l = new TH1D("hLa_dca12_qa", "Lambda QA, all masses;DCA12 (cm)", 200, 0, 0.2);
  TH1D* rl = new TH1D("hLa_rxy_qa", "Lambda QA, all masses;R_{xy} (cm)", 200, 0, 20);
  const long nk = tk.GetNtrees() ? fv::fill(tk, "K_S0", ksQa, ksLo, nullptr, nullptr, d12k, rk) : 0;
  const long nl = tl.GetNtrees() ? fv::fill(tl, "Lambda0", laQa, laLo, lbQa, lbLo, d12l, rl) : 0;
  TParameter<Long64_t>("nEventsProcessed", te.GetEntries()).Write();
  out.Write();
  out.Close();
  std::cout << "Ks QA candidates " << nk << ", Lambda QA candidates " << nl << ", events " << te.GetEntries() << std::endl;
  std::cout << "All done" << std::endl;
}
