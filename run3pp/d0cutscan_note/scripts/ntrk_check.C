// ntrk_check.C : per-candidate vs per-vertex distribution of nTracksOfVertex / nTracksOfBC (full pT>1.5 skim)
// usage (note root): root -l -b -q 'scripts/ntrk_check.C+'  -> logs/ntrk_check.txt
#include <ROOT/RDataFrame.hxx>
#include <TH1D.h>

#include <algorithm>
#include <array>
#include <cstdio>
#include <vector>

void ntrk_check()
{
  FILE *o = fopen("logs/ntrk_check.txt", "w");
  ROOT::RDataFrame d("DecayTree", "../ntuples/d0tight_full/D0tight_pT1p5_run*.root");
  auto run = d.Take<int>("runNumber");
  auto ev = d.Take<int>("eventNumber");
  auto pv = d.Take<int>("primary_vertex_ID");
  auto nv = d.Take<int>("nTracksOfVertex");
  auto nb = d.Take<int>("nTracksOfBC");
  auto npv = d.Take<int>("nPrimaryVerticesOfBC");
  auto nsi = d.Take<int>("nSiSeedMultiplicity");
  auto ndof = d.Take<unsigned int>("primary_vertex_nDoF");
  size_t n = run->size();
  std::vector<std::array<long, 7>> k(n);
  for (size_t i = 0; i < n; ++i) k[i] = {(*run)[i], (*ev)[i], (*pv)[i], (*nv)[i], (*nb)[i], (*npv)[i], (*nsi)[i]};
  std::sort(k.begin(), k.end());
  TH1D hc("hc", "", 200, 0, 200), hv("hv", "", 200, 0, 200), hb("hb", "", 400, 0, 400), hs("hs", "", 2000, 0, 20000), hp("hp", "", 20, 0, 20);
  long nVtx = 0, inconsistent = 0;
  for (size_t i = 0; i < n;)
  {
    size_t j = i;
    while (j < n && k[j][0] == k[i][0] && k[j][1] == k[i][1] && k[j][2] == k[i][2])
    {
      if (k[j][3] != k[i][3] || k[j][4] != k[i][4]) ++inconsistent;
      ++j;
    }
    hv.Fill(k[i][3]);
    hb.Fill(k[i][4]);
    hs.Fill(k[i][6]);
    hp.Fill(k[i][5]);
    ++nVtx;
    i = j;
  }
  for (size_t i = 0; i < n; ++i) hc.Fill((*nv)[i]);
  // nDoF vs nTracks: for an adaptive / simple vertex fit nDoF = 2 N - 3
  long ndofMatch = 0;
  for (size_t i = 0; i < n; ++i)
    if ((long) (*ndof)[i] == 2L * (*nv)[i] - 3) ++ndofMatch;
  auto q = [](TH1D &h, double p) { double x; h.GetQuantiles(1, &x, &p); return x; };
  fprintf(o, "candidates %zu, distinct (run, event, PV) %ld, inconsistent N within a vertex %ld\n", n, nVtx, inconsistent);
  fprintf(o, "nDoF == 2 N_vtx - 3 for %.1f%% of candidates\n", 100. * ndofMatch / n);
  for (auto *h : {&hc, &hv})
    fprintf(o, "%s N_vtx: mean %.2f, median %.0f, 90%% %.0f, frac >=20 %.2f%%, frac <=3 %.2f%%\n", h == &hc ? "per candidate" : "per vertex   ", h->GetMean(),
            q(*h, 0.5), q(*h, 0.9), 100. * h->Integral(21, 201) / h->Integral(0, 201), 100. * h->Integral(0, 4) / h->Integral(0, 201));
  fprintf(o, "per vertex N_BC: mean %.2f median %.0f; nPrimaryVerticesOfBC mean %.2f; nSiSeedMultiplicity mean %.0f median %.0f\n", hb.GetMean(), q(hb, 0.5),
          hp.GetMean(), hs.GetMean(), q(hs, 0.5));
  fprintf(o, "per vertex N_vtx histogram (N: vertices):\n");
  for (int b = 1; b <= 60; ++b) fprintf(o, "%d:%.0f%s", b - 1, hv.GetBinContent(b), b % 15 ? " " : "\n");
  fclose(o);
}
