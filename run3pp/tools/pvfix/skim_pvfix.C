// skim_pvfix.C : compact D0 skim of the d0-pvfix ntuples for the three PV versions (original, refitAll, refit)
//   root -l -b -q 'skim_pvfix.C+("<input glob>", "<output file>")'
// Keeps a candidate if both daughters have pT > 0.5 and it passes |DCA_PV| > 60 um (both daughters) and DIRA > 0.95
// with respect to ANY of the three PVs. Stores the PV-dependent quantities of all three versions.
#include <ROOT/RDataFrame.hxx>

#include <cstdio>
#include <string>
#include <vector>

void skim_pvfix(const char *in, const char *out)
{
  ROOT::EnableImplicitMT(4);
  ROOT::RDataFrame d("DecayTree", in);
  auto pass = [](const std::string &s) {
    return "(abs(track_1_PV_DCA" + s + ")>0.006 && abs(track_2_PV_DCA" + s + ")>0.006 && D0_DIRA" + s + ">0.95)";
  };
  std::string cut = "track_1_pT>0.5 && track_2_pT>0.5 && (" + pass("") + " || (primary_vertex_refitAll_status==1 && " + pass("_PVrefitAll") +
                    ") || (primary_vertex_refit_status==1 && " + pass("_PVrefit") + "))";
  std::vector<std::string> cols = {"runNumber", "eventNumber", "D0_mass", "D0_massErr", "D0_pT", "D0_pseudorapidity", "D0_phi", "D0_chi2", "D0_nDoF",
                                   "track_1_track_2_DCA", "nTracksOfVertex", "nTracksOfVertex_noDaughters", "nPrimaryVerticesOfBC", "primary_vertex_ID",
                                   "primary_vertex_crossing", "primary_vertex_refit_status", "primary_vertex_refitAll_status", "primary_vertex_refit_nTracks",
                                   "primary_vertex_refitAll_nTracks"};
  for (int i = 1; i <= 2; ++i)
    for (const char *v : {"pT", "pseudorapidity", "phi", "charge", "PDG_ID", "track_ID", "calculated_dEdx", "inPV"}) cols.push_back(Form("track_%d_%s", i, v));
  for (const std::string s : {"", "_PVrefitAll", "_PVrefit"})
  {
    for (const char *v : {"DIRA", "decayLength", "decayLengthErr", "PV_DCA"}) cols.push_back(std::string("D0_") + v + s);
    for (int i = 1; i <= 2; ++i)
      for (const char *v : {"PV_DCA", "PV_DCA_xy", "PV_DCA_sig"}) cols.push_back(Form("track_%d_%s%s", i, v, s.c_str()));
  }
  // primary vertex positions of the three versions (for vertex-shift studies)
  for (const char *v : {"primary_vertex_x", "primary_vertex_y", "primary_vertex_z", "primary_vertex_refitAll_x", "primary_vertex_refitAll_y",
                        "primary_vertex_refitAll_z", "primary_vertex_refit_x", "primary_vertex_refit_y", "primary_vertex_refit_z"})
    cols.push_back(v);
  ROOT::RDF::RSnapshotOptions o;
  o.fCompressionAlgorithm = ROOT::kZSTD;
  o.fCompressionLevel = 5;
  auto n = d.Count();
  auto s = d.Filter(cut).Snapshot("DecayTree", out, cols, o);
  auto ns = (*s).Count();
  printf("skim_pvfix: %llu -> %llu candidates (%zu branches)\n", *n, *ns, cols.size());
}
