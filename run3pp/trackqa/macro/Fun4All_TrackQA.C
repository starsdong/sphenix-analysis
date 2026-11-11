/*
 * Event/track QA ntuple from run3pp DST_TRKR_TRACKS (ana581).
 *
 * Reads the DST (no geometry, no CDB, no reconstruction) and writes with TrackQATree:
 *   <outdir>/TrackQA_<run>_<seg>.root   trees "event" (all events) and "track" (all tracks)
 * Optional (trimDST=true): a trimmed DST with only the nodes needed for track-level analysis
 *   <outdir>/DST_TRKR_TRACKS_TRIM-<run>-<seg>.root
 *   kept: Sync, EventHeader, GL1RAWHIT, SvtxTrackMap, SvtxVertexMap
 *   dropped: TRKR_CLUSTER_SEED (clusters), Silicon/Tpc/SvtxTrackSeedContainer (seeds)
 *
 * Usage: root.exe -b -q 'Fun4All_TrackQA.C(nEvents, "input.root or list.txt", "outdir", trimDST, minPt)'
 *   nEvents <= 0 : all events
 */

#include <trackqatree/TrackQATree.h>

#include <fun4all/Fun4AllDstInputManager.h>
#include <fun4all/Fun4AllDstOutputManager.h>
#include <fun4all/Fun4AllServer.h>
#include <fun4all/Fun4AllUtils.h>
#include <phool/recoConsts.h>

#include <TSystem.h>

#include <fstream>
#include <iomanip>
#include <sstream>

R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libffarawobjects.so)
R__LOAD_LIBRARY(libtrack_io.so)
R__LOAD_LIBRARY(libtrackbase_historic_io.so)
R__LOAD_LIBRARY(libglobalvertex_io.so)
R__LOAD_LIBRARY(libtrackqatree.so)
// dictionaries for the geometry objects stored in the DST run tree (otherwise ROOT read errors
// and broken run nodes in the trimmed DST)
R__LOAD_LIBRARY(libg4detectors_io.so)
R__LOAD_LIBRARY(libmvtx_io.so)
R__LOAD_LIBRARY(libintt_io.so)
R__LOAD_LIBRARY(libtpc_io.so)
R__LOAD_LIBRARY(libmicromegas_io.so)

void Fun4All_TrackQA(const int nEvents = 0,
                     const std::string &input = "input.root",
                     const std::string &outdir = "./",
                     const bool trimDST = false,
                     const float minPt = 0.)
{
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
  std::stringstream tag;
  tag << std::setw(8) << std::setfill('0') << runseg.first << "_" << std::setw(5) << std::setfill('0') << runseg.second;
  gSystem->mkdir(outdir.c_str(), true);
  const std::string qaout = outdir + "/TrackQA_" + tag.str() + ".root";
  // trimmed DST keeps the production "-<run>-<seg>" pattern so Fun4AllUtils::GetRunSegment works on it
  std::stringstream dtag;
  dtag << std::setw(8) << std::setfill('0') << runseg.first << "-" << std::setw(5) << std::setfill('0') << runseg.second;
  const std::string dstout = outdir + "/DST_TRKR_TRACKS_TRIM-" + dtag.str() + ".root";

  auto *se = Fun4AllServer::instance();
  se->Verbosity(1);
  recoConsts::instance()->set_IntFlag("RUNNUMBER", runseg.first);

  auto *in = new Fun4AllDstInputManager("TrackInputManager");
  for (const auto &f : files)
  {
    std::cout << "Adding input " << f << std::endl;
    in->AddFile(f);
  }
  se->registerInputManager(in);

  auto *qa = new TrackQATree("TrackQATree", qaout);
  qa->setMinPt(minPt);
  se->registerSubsystem(qa);

  if (trimDST)
  {
    auto *out = new Fun4AllDstOutputManager("TrimOut", dstout);
    out->AddNode("Sync");
    out->AddNode("EventHeader");
    out->AddNode("GL1RAWHIT");
    out->AddNode("SvtxTrackMap");
    out->AddNode("SvtxVertexMap");
    se->registerOutputManager(out);
  }

  se->run(nEvents > 0 ? nEvents : 0);
  se->End();
  se->PrintTimer();

  std::vector<std::string> outputs = {qaout};
  if (trimDST) outputs.push_back(dstout);
  int nmissing = 0;
  for (const auto &o : outputs)
  {
    bool ok = !gSystem->AccessPathName(o.c_str());
    std::cout << (ok ? "OUTPUT " : "MISSING ") << o << std::endl;
    if (!ok) ++nmissing;
  }
  delete se;
  if (nmissing == 0) std::cout << "All done" << std::endl;
  gSystem->Exit(nmissing == 0 ? 0 : 1);
}
