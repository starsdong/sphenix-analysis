/*
 * run3pp KFParticle V0 / open-charm reconstruction on DST_TRKR_TRACKS
 *   K_S0   -> pi+ pi-
 *   Lambda -> p pi-   (+ c.c.)
 *   D0     -> K- pi+  (+ c.c.)
 *
 * Loose topological cuts; every topological variable KFParticle provides is written to
 * one flat ntuple (tree "DecayTree") per channel, plus a per-event companion tree
 * (EventInfoTree, tree "EventInfo": all events, vertex list with beam crossings).
 *
 * Crossing handling (streaming readout, many crossings per event):
 *   default channels   : daughters must share one crossing AND the PV must have that crossing
 *                        (requireBunchCrossingMatch + requireTrackVertexBunchCrossingMatch)
 *   crossingStudy=true : extra Ks/Lambda instances ("*_noXmatch") with both crossing
 *                        requirements switched off, to check what the crossing choice does.
 *
 * Based on macros/TrackingProduction/Fun4All_KShortReco.C (macros commit 87f5cb5e = ana.581).
 *
 * Ntuple notes:
 *   - primary_vertex_ID holds the BEAM CROSSING of the chosen PV (KFParticle getPVID), not the map key.
 *   - track_N_bunch_crossing is relative to the event's GL1 BCO; TPC-only tracks always have 0.
 *   - Collision_BCO = GL1_BCO + track_1 crossing: the same collision can be reconstructed in
 *     consecutive (overlapping streaming) events -> de-duplicate offline with Collision_BCO.
 * Cut sets (last argument):
 *   "loose" (default) : very loose pre-cuts, for small-sample diagnostics
 *   "prod"            : trimmed for mass production:
 *                       Ks, Lambda : daughter 3D DCA to PV > 0.05 cm, DIRA > 0.9
 *                       D0         : daughter 3D DCA to PV > 0.003 cm (30 um), DIRA > 0.9
 *                       (use with crossingStudy=false: no *_noXmatch output)
 *                       output format slimmed: no trigger block, redundant momentum and all covariance branches removed,
 *                       ZSTD recompression (tools/slim_kfp.C)
 *   KFParticle's track DCA cut uses the minimum 3D DCA over the PVs of the track's own crossing
 *   (all PVs for TPC-only tracks).
 * Usage: root.exe -b -q 'Fun4All_run3pp_KFP.C(nEvents, "input.root or list.txt", "outdir", crossingStudy, nSkip, "Ks,Lambda,D0", "loose|prod")'
 */

#include <GlobalVariables.C>

#include <G4_ActsGeom.C>
#include <G4_Magnet.C>
#include <Trkr_RecoInit.C>
#include <Trkr_TpcReadoutInit.C>

#include <kfparticle_sphenix/KFParticle_sPHENIX.h>
#include <eventinfotree/EventInfoTree.h>
#include <calotrigger/TriggerRunInfoReco.h>

#include <ffamodules/CDBInterface.h>
#include <fun4all/Fun4AllDstInputManager.h>
#include <fun4all/Fun4AllRunNodeInputManager.h>
#include <fun4all/Fun4AllServer.h>
#include <fun4all/Fun4AllUtils.h>
#include <phool/recoConsts.h>

#include <TSystem.h>
#include <TFile.h>
#include <TParameter.h>

#include "slim_kfp.C"  // copied next to this macro by run_kfp.sh (from tools/)

#include <cfloat>
#include <fstream>
#include <iomanip>
#include <sstream>

R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libffamodules.so)
R__LOAD_LIBRARY(libphool.so)
R__LOAD_LIBRARY(libcdbobjects.so)
R__LOAD_LIBRARY(libcalotrigger.so)
R__LOAD_LIBRARY(libkfparticle_sphenix.so)
R__LOAD_LIBRARY(libeventinfotree.so)

namespace KFPrun3pp
{
  bool getTriggerInfo = true;   // 64 trigger-bit branches per candidate (needs TriggerRunInfo)
  bool getDetectorInfo = true;  // per-daughter MVTX/INTT/TPC/TPOT nHits / nStates
  bool getDetailedTracking = false;  // per-cluster vectors (layer, residuals, stave/chip/sector IDs): large
  bool getdEdx = true;          // per-daughter TPC dE/dx
  bool getAllPVInfo = false;    // vectors with DCA to every PV: large, off by default
  bool prod = false;            // production cut set (see header)
  int verbosity = 0;
}  // namespace KFPrun3pp

// settings shared by all channels
void kfpCommon(KFParticle_sPHENIX *kfp, const std::string &outfile, bool xMatch)
{
  kfp->Verbosity(KFPrun3pp::verbosity);
  kfp->usePID(false);  // loose: no dE/dx PID cut, dE/dx is stored instead
  kfp->dontUseGlobalVertex(true);
  kfp->requireBunchCrossingMatch(xMatch);             // all daughters from the same crossing
  kfp->requireTrackVertexBunchCrossingMatch(xMatch);  // PV crossing == daughter crossing
  kfp->getAllPVInfo(KFPrun3pp::getAllPVInfo);
  kfp->allowZeroMassTracks();
  kfp->use2Dmatching(false);
  kfp->getTriggerInfo(KFPrun3pp::getTriggerInfo);
  kfp->getDetectorInfo(KFPrun3pp::getDetectorInfo);
  kfp->GetDetailedTracking(KFPrun3pp::getDetailedTracking);  // default in KFP is true!
  kfp->get_dEdx_info(KFPrun3pp::getdEdx);
  kfp->saveDST(false);
  kfp->saveParticleContainer(false);
  kfp->magFieldFile("FIELDMAP_TRACKING");
  kfp->constrainToPrimaryVertex(true);  // pick a PV, compute pointing / decay-length w.r.t. it

  // remove all lifetime / significance requirements (stored, cut offline)
  kfp->setFlightDistancechi2(-1.);
  kfp->setDecayLengthRange(-FLT_MAX, FLT_MAX);
  kfp->setDecayLengthRange_XY(-FLT_MAX, FLT_MAX);
  kfp->setDecayTimeRange(-FLT_MAX, FLT_MAX);
  kfp->setDecayTimeRange_XY(-FLT_MAX, FLT_MAX);
  kfp->setMinDecayTimeSignificance(-FLT_MAX);
  kfp->setMinDecayLengthSignificance(-FLT_MAX);
  kfp->setMinDecayLengthSignificance_XY(-FLT_MAX);
  kfp->setMotherPV_DCA_StdDev(FLT_MAX);
  kfp->setMaximumDaughterDCA_XY(FLT_MAX);
  kfp->setMinimumTrackPV_DCA(-1.);
  kfp->setMinimumTrackPV_DCA_StdDev(-1.);
  kfp->setMaximumTrackchi2nDOF(FLT_MAX);
  kfp->setMotherPT(0);

  kfp->setOutputName(outfile);
}

KFParticle_sPHENIX *makeKshort(const std::string &name, const std::string &outfile, bool xMatch)
{
  auto *kfp = new KFParticle_sPHENIX(name);
  kfp->setDecayDescriptor("K_S0 -> pi^+ pi^-");
  kfpCommon(kfp, outfile, xMatch);
  kfp->setMinimumTrackPT(0.1);
  kfp->setMinMVTXhits(0);
  kfp->setMinINTThits(0);
  kfp->setMinTPChits(20);
  kfp->setMaximumDaughterDCA(0.5);  // cm
  kfp->setMaximumVertexchi2nDOF(FLT_MAX);  // [loosecov] no vertex-chi2 cut
  kfp->setMinDIRA(0.8);
  kfp->setMinimumMass(0.40);
  kfp->setMaximumMass(0.60);
  kfp->setMaximumMotherVertexVolume(FLT_MAX);  // [loosecov] no vertex-volume cut
  if (KFPrun3pp::prod)
  {
    kfp->setMinimumTrackPV_DCA(0.05);  // cm
    kfp->setMinDIRA(0.9);
  }
  return kfp;
}

KFParticle_sPHENIX *makeLambda(const std::string &name, const std::string &outfile, bool xMatch)
{
  auto *kfp = new KFParticle_sPHENIX(name);
  kfp->setDecayDescriptor("[Lambda0 -> proton^+ pi^-]cc");
  kfpCommon(kfp, outfile, xMatch);
  kfp->setMinimumTrackPT(0.1);
  kfp->setMinMVTXhits(0);
  kfp->setMinINTThits(0);
  kfp->setMinTPChits(20);
  kfp->setMaximumDaughterDCA(0.5);  // cm
  kfp->setMaximumVertexchi2nDOF(FLT_MAX);  // [loosecov] no vertex-chi2 cut
  kfp->setMinDIRA(0.8);
  kfp->setMinimumMass(1.05);
  kfp->setMaximumMass(1.20);
  kfp->setMaximumMotherVertexVolume(FLT_MAX);  // [loosecov] no vertex-volume cut
  if (KFPrun3pp::prod)
  {
    kfp->setMinimumTrackPV_DCA(0.05);  // cm
    kfp->setMinDIRA(0.9);
  }
  return kfp;
}

KFParticle_sPHENIX *makeD0(const std::string &name, const std::string &outfile, bool xMatch)
{
  // NOTE: with "[...]cc" KFParticle keeps one mass hypothesis per track pair (best PV-DCA chi2);
  // the swapped (K<->pi) mass can be recomputed offline from the stored daughter momenta.
  auto *kfp = new KFParticle_sPHENIX(name);
  kfp->setDecayDescriptor("[D0 -> K^- pi^+]cc");
  kfpCommon(kfp, outfile, xMatch);
  kfp->setMinimumTrackPT(0.2);
  kfp->setMinMVTXhits(2);  // need MVTX pointing for charm
  kfp->setMinINTThits(0);
  kfp->setMinTPChits(20);
  kfp->setMaximumDaughterDCA(0.05);  // 500 um
  kfp->setMaximumVertexchi2nDOF(FLT_MAX);  // [loosecov] no vertex-chi2 cut
  kfp->setMinDIRA(0.0);
  kfp->setMinimumMass(1.60);
  kfp->setMaximumMass(2.10);
  kfp->setMaximumMotherVertexVolume(FLT_MAX);  // [loosecov] no vertex-volume cut
  if (KFPrun3pp::prod)
  {
    kfp->setMinimumTrackPV_DCA(0.003);  // cm = 30 um
    kfp->setMinDIRA(0.9);
  }
  return kfp;
}

// ---------------------------------------------------------------------------------------------
// 3-prong charm hadrons: D+ -> K- pi+ pi+, Ds+ -> phi(K+K-) pi+, Lambda_c+ -> p K- pi+ (and c.c.)
// All cuts in one struct so that they are printed into the log and easy to tune.
struct CharmCuts
{
  std::string name, descriptor;
  float trkPt;        // GeV/c, every daughter
  int minMVTX;        // MVTX states per daughter
  int minTPC;         // TPC states per daughter
  float trkPVDCA;     // cm, min 3D DCA of every daughter to the PVs of its crossing
  float maxPairDCA;   // cm, max DCA between any two daughters
  float maxVtxChi2;   // secondary-vertex chi2/ndf
  float minDIRA;      // cos(pointing angle) w.r.t. the PV
  float mLo, mHi;     // GeV/c^2, mother mass window
  float maxVtxVol;    // cm^3
  float phiLo, phiHi; // GeV/c^2, intermediate phi mass window (Ds only)
};

CharmCuts charmCuts(const std::string &ch)
{
  //                 name     descriptor                                  pT   MVTX TPC  trkDCA  pairDCA chi2 DIRA  mLo   mHi   vol  phiLo phiHi
  if (ch == "Dplus") return {"Dplus", "[D+ -> K^- pi^+ pi^+]cc",            0.3, 2,  20,  0.005,  0.03,   10,  0.95, 1.70, 2.05, 0.1, 0, 0};
  if (ch == "Ds")    return {"Ds", "[D_s+ -> {phi -> K^+ K^-} pi^+]cc",      0.3, 2,  20,  0.003,  0.03,   10,  0.95, 1.75, 2.10, 0.1, 0.99, 1.05};
  if (ch == "Lc")    return {"Lc", "[Lambda_c+ -> proton^+ K^- pi^+]cc",     0.3, 2,  20,  0.002,  0.02,   10,  0.95, 2.15, 2.42, 0.1, 0, 0};
  std::cout << "charmCuts: unknown channel " << ch << std::endl;
  gSystem->Exit(1);
  return {};
}

void printCharmCuts(const CharmCuts &c)
{
  std::cout << "CHARMCUTS " << c.name << " | " << c.descriptor << " | trk pT>" << c.trkPt << " MVTX>=" << c.minMVTX
            << " TPC>=" << c.minTPC << " trkPVDCA>" << c.trkPVDCA << " | pairDCA<" << c.maxPairDCA << " vtxChi2/ndf<"
            << c.maxVtxChi2 << " DIRA>" << c.minDIRA << " mass " << c.mLo << "-" << c.mHi << " vtxVol<" << c.maxVtxVol;
  if (c.phiHi > 0) std::cout << " | phi mass " << c.phiLo << "-" << c.phiHi;
  std::cout << std::endl;
}

KFParticle_sPHENIX *makeCharm3(const std::string &ch, const std::string &outfile)
{
  const CharmCuts c = charmCuts(ch);
  printCharmCuts(c);
  auto *kfp = new KFParticle_sPHENIX(c.name + "_reco");
  kfp->setDecayDescriptor(c.descriptor);
  kfpCommon(kfp, outfile, true);
  kfp->setMinimumTrackPT(c.trkPt);
  kfp->setMinMVTXhits(c.minMVTX);
  kfp->setMinINTThits(0);
  kfp->setMinTPChits(c.minTPC);
  kfp->setMinimumTrackPV_DCA(c.trkPVDCA);
  kfp->setMaximumDaughterDCA(c.maxPairDCA);
  kfp->setMaximumVertexchi2nDOF(c.maxVtxChi2);
  kfp->setMinDIRA(c.minDIRA);
  kfp->setMinimumMass(c.mLo);
  kfp->setMaximumMass(c.mHi);
  kfp->setMaximumMotherVertexVolume(c.maxVtxVol);
  if (c.phiHi > 0)
  {
    // phi decays at the Ds vertex: no pointing / flight / PV-DCA requirement on it
    kfp->constrainIntermediateMasses(false);
    kfp->setIntermediateMassRange({{c.phiLo, c.phiHi}});
    kfp->setIntermediateMinPT({0.});
    kfp->setIntermediateMinDIRA({-1.1});
    kfp->setIntermediateMinFDchi2({-1.});
    kfp->setIntermediateMaxVertexVolume({c.maxVtxVol});
    kfp->setIntermediatePV_DCARange({{-1., FLT_MAX}});
    kfp->setIntermediatePV_DCARange_XY({{-1., FLT_MAX}});
    kfp->setIntermediatePV_DCA_StdDevRange({{-1., FLT_MAX}});
    kfp->setIntermediatePV_DCA_StdDevRange_XY({{-1., FLT_MAX}});
  }
  return kfp;
}

void Fun4All_run3pp_KFP(const int nEvents = 100,
                        const std::string &input = "input.root",
                        const std::string &outdir = "./",
                        const bool crossingStudy = false,
                        const int nSkip = 0,
                        const std::string &channels = "Ks,Lambda,D0",
                        const std::string &cutset = "loose")
{
  if (cutset != "loose" && cutset != "prod")
  {
    std::cout << "Unknown cut set " << cutset << std::endl;
    gSystem->Exit(1);
  }
  KFPrun3pp::prod = (cutset == "prod");
  if (KFPrun3pp::prod)
  {
    KFPrun3pp::getTriggerInfo = false;  // trigger bits per event are in EventInfo (GL1_scaled/GL1_live)
  }
  std::cout << "KFP cut set: " << cutset << std::endl;
  // ---- inputs: a single DST or a text list of DSTs (first file defines run/segment)
  std::vector<std::string> files;
  if (input.size() > 5 && input.substr(input.size() - 5) == ".root")
  {
    files.push_back(input);
  }
  else
  {
    std::ifstream ifs(input);
    std::string line;
    while (std::getline(ifs, line))
    {
      if (!line.empty()) files.push_back(line);
    }
  }
  if (files.empty())
  {
    std::cout << "No input files" << std::endl;
    gSystem->Exit(1);
  }

  std::pair<int, int> runseg = Fun4AllUtils::GetRunSegment(files[0]);
  const int runnumber = runseg.first;
  const int segment = runseg.second;
  std::stringstream tag;
  tag << std::setw(8) << std::setfill('0') << runnumber << "_" << std::setw(5) << std::setfill('0') << segment
      << "_" << std::setw(5) << std::setfill('0') << nSkip;

  gSystem->mkdir(outdir.c_str(), true);
  auto outname = [&](const std::string &ch) { return outdir + "/KFP_" + ch + "_" + tag.str() + ".root"; };

  auto *rc = recoConsts::instance();
  rc->set_IntFlag("RUNNUMBER", runnumber);
  rc->set_IntFlag("RUNSEGMENT", segment);
  Enable::CDB = true;
  rc->set_StringFlag("CDB_GLOBALTAG", "newcdbtag");
  rc->set_uint64Flag("TIMESTAMP", runnumber);

  // geometry / field / TPC readout config as in the production KShortReco macro
  std::string geofile = CDBInterface::instance()->getUrl("Tracking_Geometry");
  TpcReadoutInit(runnumber);
  G4TRACKING::SC_CALIBMODE = false;
  Enable::MVTX_APPLYMISALIGNMENT = true;
  ACTSGEOM::mvtx_applymisalignment = Enable::MVTX_APPLYMISALIGNMENT;
  TRACKING::streaming_mode = true;

  auto *se = Fun4AllServer::instance();
  se->Verbosity(1);

  auto *ingeo = new Fun4AllRunNodeInputManager("GeoIn");
  ingeo->AddFile(geofile);
  se->registerInputManager(ingeo);

  TrackingInit();

  auto *in = new Fun4AllDstInputManager("TrackInputManager");
  for (const auto &f : files)
  {
    std::cout << "Adding input " << f << std::endl;
    in->AddFile(f);
  }
  se->registerInputManager(in);

  if (KFPrun3pp::getTriggerInfo)
  {
    se->registerSubsystem(new TriggerRunInfoReco());
  }

  auto want = [&](const std::string &ch) { return ("," + channels + ",").find("," + ch + ",") != std::string::npos; };
  // channel token "noEventInfo": no per-event EventInfo tree (space saving); the number of processed events
  // is then only stored as TParameter<Long64_t> "nEventsProcessed" in each output file (see below)
  const bool evInfo = !want("noEventInfo");
  std::vector<std::string> outputs;
  if (evInfo)
  {
    se->registerSubsystem(new EventInfoTree("EventInfoTree", outname("EventInfo")));
    outputs.push_back(outname("EventInfo"));
  }
  if (want("Ks"))
  {
    se->registerSubsystem(makeKshort("Ks_reco", outname("Ks"), true));
    outputs.push_back(outname("Ks"));
  }
  if (want("Lambda"))
  {
    se->registerSubsystem(makeLambda("Lambda_reco", outname("Lambda"), true));
    outputs.push_back(outname("Lambda"));
  }
  if (want("D0"))
  {
    se->registerSubsystem(makeD0("D0_reco", outname("D0"), true));
    outputs.push_back(outname("D0"));
  }
  for (const char *ch3 : {"Dplus", "Ds", "Lc"})
  {
    if (want(ch3))
    {
      se->registerSubsystem(makeCharm3(ch3, outname(ch3)));
      outputs.push_back(outname(ch3));
    }
  }
  if (crossingStudy)
  {
    se->registerSubsystem(makeKshort("Ks_noXmatch", outname("Ks_noXmatch"), false));
    se->registerSubsystem(makeLambda("Lambda_noXmatch", outname("Lambda_noXmatch"), false));
    outputs.push_back(outname("Ks_noXmatch"));
    outputs.push_back(outname("Lambda_noXmatch"));
  }

  se->skip(nSkip);
  se->run(nEvents);
  se->End();
  se->PrintTimer();

  if (KFPrun3pp::prod)
  {
    for (const auto &o : outputs)
    {
      if (o.find("EventInfo") == std::string::npos && !gSystem->AccessPathName(o.c_str()) && !slimDecayTree(o))
      {
        std::cout << "slimming FAILED for " << o << std::endl;
        gSystem->Exit(1);
      }
    }
  }

  // number of processed events, stored in every candidate file (hadd sums TParameter<Long64_t>)
  const Long64_t nProcessed = se->EventCounter();
  std::cout << "nEventsProcessed " << nProcessed << std::endl;
  for (const auto &o : outputs)
  {
    if (o.find("EventInfo") != std::string::npos || gSystem->AccessPathName(o.c_str())) continue;
    TFile fu(o.c_str(), "UPDATE");
    TParameter<Long64_t> np("nEventsProcessed", nProcessed);
    np.Write();
    fu.Close();
  }

  int nmissing = 0;
  for (const auto &o : outputs)
  {
    bool ok = !gSystem->AccessPathName(o.c_str());
    std::cout << (ok ? "OUTPUT " : "MISSING ") << o << std::endl;
    if (!ok) ++nmissing;
  }
  delete se;
  if (nmissing == 0)
  {
    std::cout << "All done" << std::endl;
  }
  gSystem->Exit(nmissing == 0 ? 0 : 1);
}
