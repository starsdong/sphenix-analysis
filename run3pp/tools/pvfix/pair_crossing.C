// pair_crossing.C : what track pairs enter the D0 reconstruction? Crossing match of the daughters, PV crossing,
// crossing values, INTT content, number of PVs and tracks in the crossing, daughter attachment to PVs.
//   root -l -b -q 'pair_crossing.C+("<file glob>")'
#include <ROOT/RDataFrame.hxx>
#include <TH1D.h>

#include <cstdio>

void pair_crossing(const char *in)
{
  ROOT::RDataFrame d("DecayTree", in);
  auto n = d.Count();
  auto same = d.Filter("track_1_bunch_crossing==track_2_bunch_crossing").Count();
  auto pvx = d.Filter("track_1_bunch_crossing==primary_vertex_crossing").Count();
  auto x0 = d.Filter("track_1_bunch_crossing==0").Count();
  auto neg = d.Filter("track_1_bunch_crossing<0").Count();
  auto big = d.Filter("track_1_bunch_crossing>400").Count();
  auto bothI = d.Filter("track_1_INTT_nHits>0 && track_2_INTT_nHits>0").Count();
  auto noneI = d.Filter("track_1_INTT_nHits==0 && track_2_INTT_nHits==0").Count();
  auto oneI = d.Filter("(track_1_INTT_nHits>0) != (track_2_INTT_nHits>0)").Count();
  auto both = d.Filter("track_1_inPV==1 && track_2_inPV==1").Count();
  auto none = d.Filter("track_1_inPV==0 && track_2_inPV==0").Count();
  auto multiPV = d.Filter("nPrimaryVerticesOfBC>1").Count();
  auto hPV = d.Histo1D({"pv", "", 10, 0, 10}, "nPrimaryVerticesOfBC");
  auto hT = d.Histo1D({"t", "", 400, 0, 400}, "nTracksOfBC");
  auto hV = d.Histo1D({"v", "", 400, 0, 400}, "nTracksOfVertex");
  auto hx = d.Histo1D({"x", "", 700, -100, 600}, "track_1_bunch_crossing");
  double N = *n;
  printf("candidates %llu\n", *n);
  printf("daughters in the same crossing: %.4f | daughter crossing == PV crossing: %.4f\n", *same / N, *pvx / N);
  printf("crossing == 0: %.3f, < 0: %.3f, > 400: %.3f\n", *x0 / N, *neg / N, *big / N);
  double q[5], p[5] = {0.05, 0.25, 0.5, 0.75, 0.95};
  hx->GetQuantiles(5, q, p);
  printf("crossing quantiles 5/25/50/75/95%%: %.0f %.0f %.0f %.0f %.0f\n", q[0], q[1], q[2], q[3], q[4]);
  printf("INTT hits: both daughters %.3f, one %.3f, neither %.3f\n", *bothI / N, *oneI / N, *noneI / N);
  printf("daughters attached to the candidate's PV: both %.3f, neither %.3f\n", *both / N, *none / N);
  printf("PVs in the crossing:");
  for (int i = 1; i <= 6; ++i) printf(" %d:%.3f", i - 1, hPV->GetBinContent(i) / N);
  printf("  (>1: %.3f)\n", *multiPV / N);
  hT->GetQuantiles(5, q, p);
  printf("tracks in the crossing (nTracksOfBC) quantiles 5/25/50/75/95%%: %.0f %.0f %.0f %.0f %.0f, mean %.1f\n", q[0], q[1], q[2], q[3], q[4], hT->GetMean());
  hV->GetQuantiles(5, q, p);
  printf("tracks of the PV (nTracksOfVertex) quantiles: %.0f %.0f %.0f %.0f %.0f, mean %.1f\n", q[0], q[1], q[2], q[3], q[4], hV->GetMean());
}
