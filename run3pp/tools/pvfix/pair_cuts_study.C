// pair_cuts_study.C : effect of the pair-quality cuts on the v2 D0 skim (original PV, L < 1 mm, both |DCA| > 60 um).
// For the pointing (DIRA > 0.995, L/sigma > 2) and control (DIRA 0.95-0.98) samples: pass fractions of the INTT and
// same-collision cuts, split by crossings with one or several PVs, and the |dz| distribution of non-attached daughters.
//   root -l -b -q 'pair_cuts_study.C+("<skim>")'
#include <ROOT/RDataFrame.hxx>
#include <TH1D.h>

#include <cstdio>
#include <string>

void pair_cuts_study(const char *in)
{
  ROOT::EnableImplicitMT(6);
  ROOT::RDataFrame d0("DecayTree", in);
  auto d = d0.Filter("pvfix_passOriginal==1 && abs(track_1_PV_DCA)>0.006 && abs(track_2_PV_DCA)>0.006 && D0_decayLength<0.1")
               .Define("dz1", "float a=track_1_PV_DCA, b=track_1_PV_DCA_xy; return (float) std::sqrt(std::max(0.f, a*a-b*b));")
               .Define("dz2", "float a=track_2_PV_DCA, b=track_2_PV_DCA_xy; return (float) std::sqrt(std::max(0.f, a*a-b*b));")
               .Define("intt", "track_1_INTT_nHits>0 && track_2_INTT_nHits>0")
               .Define("coll", "(track_1_inPV==1 || dz1<0.1) && (track_2_inPV==1 || dz2<0.1)")
               .Define("dzNA", "ROOT::RVecF v; if (track_1_inPV!=1) v.push_back(dz1); if (track_2_inPV!=1) v.push_back(dz2); return v;");
  const char *selName[2] = {"pointing (DIRA>0.995, L/s>2)", "control (DIRA 0.95-0.98)"};
  const char *sel[2] = {"D0_DIRA>0.995 && D0_decayLength/D0_decayLengthErr>2", "D0_DIRA>0.95 && D0_DIRA<0.98"};
  for (int s = 0; s < 2; ++s)
  {
    printf("== %s\n", selName[s]);
    for (const char *pv : {"nPrimaryVerticesOfBC>=1", "nPrimaryVerticesOfBC==1", "nPrimaryVerticesOfBC>1"})
    {
      auto x = d.Filter(sel[s]).Filter(pv);
      auto n = x.Count(), ni = x.Filter("intt").Count(), nc = x.Filter("coll").Count(), nb = x.Filter("intt && coll").Count();
      auto none = x.Filter("track_1_inPV!=1 && track_2_inPV!=1").Count(), one = x.Filter("(track_1_inPV==1) != (track_2_inPV==1)").Count();
      auto h = x.Histo1D({"h", "", 2000, 0, 20}, "dzNA");
      double N = *n, q[3], p[3] = {0.5, 0.75, 0.9};
      h->GetQuantiles(3, q, p);
      printf("  %-24s N %9.0f | INTT %.3f | same-coll %.3f | both %.3f | daughters attached: one %.3f, none %.3f | non-attached |dz|: "
             "median %.3f, 75%% %.3f, 90%% %.3f cm, > 1 mm %.3f\n",
             pv, N, *ni / N, *nc / N, *nb / N, *one / N, *none / N, q[0], q[1], q[2], h->Integral(h->FindFixBin(0.1), h->GetNbinsX() + 1) / std::max(1.0, h->GetEntries()));
    }
  }
}
