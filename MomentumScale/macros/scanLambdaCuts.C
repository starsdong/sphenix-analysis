// scanLambdaCuts.C -- S/B scan of tighter Lambda topological cuts on the KFParticle Lambda ntuple (run 79509 seg 0-30).
//   root -l -b -q 'macros/scanLambdaCuts.C+("data/loosecov/KFP_Lambda_run79509_seg000-030_loosecov.root")'
// Base = run3pp QA selection: |PVz|<10, |eta_d|<1.1, DIRA>0.995, Lxy>0.5, DCA12<0.2, |DCAxy(pi)|,|DCAxy(p)|>0.02 cm.
// Only cuts DDSim can model (DIRA, L_xy, daughter |DCA_xy|) are scanned. Mass = KFP mass (scan only).
// Signal: counts in 1.1097-1.1217 (+-6 MeV) minus linear background from side bands 1.095-1.104, 1.128-1.137 (scaled by width).
// Prints SCAN rows: piDca pDca lxy dira | per pT bin S B ; integrated S B S/B Seff=S^2/(S+B)
#include <TFile.h>
#include <TTree.h>
#include <cmath>
#include <cstdio>
#include <vector>

void scanLambdaCuts(const char* file)
{
  TFile f(file); TTree* t = (TTree*)f.Get("DecayTree");
  float m, pt, dira, lxy, d12, pvz, e1, e2, x1, x2;
  t->SetBranchStatus("*", 0);
  const char* br[] = {"Lambda0_mass", "Lambda0_pT", "Lambda0_DIRA", "Lambda0_decayLength_xy", "track_1_track_2_DCA",
                      "primary_vertex_z", "track_1_pseudorapidity", "track_2_pseudorapidity", "track_1_PV_DCA_xy", "track_2_PV_DCA_xy"};
  float* ad[] = {&m, &pt, &dira, &lxy, &d12, &pvz, &e1, &e2, &x1, &x2};
  for (int i = 0; i < 10; ++i) { t->SetBranchStatus(br[i], 1); t->SetBranchAddress(br[i], ad[i]); }
  struct C { float m, pt, dira, lxy, dpi, dp; };
  std::vector<C> v;
  for (Long64_t i = 0; i < t->GetEntries(); ++i) {
    t->GetEntry(i);
    if (std::fabs(pvz) >= 10 || std::fabs(e1) >= 1.1 || std::fabs(e2) >= 1.1) continue;
    if (!(dira > 0.995 && lxy > 0.5 && d12 < 0.2 && std::fabs(x1) > 0.02 && std::fabs(x2) > 0.02)) continue;
    if (m < 1.09 || m > 1.14) continue;
    v.push_back({m, pt, dira, lxy, std::fabs(x1), std::fabs(x2)});   // track_1 = pion, track_2 = proton
  }
  printf("base candidates in 1.09-1.14: %zu\n", v.size());
  const std::vector<double> ptb = {0.5, 0.9, 1.2, 1.6, 2.4, 4.0};
  const double sLo = 1.1097, sHi = 1.1217, b1 = 1.095, b2 = 1.104, b3 = 1.128, b4 = 1.137;
  const double wS = sHi - sLo, wB = (b2 - b1) + (b4 - b3);
  printf("SCANHDR piDca pDca lxy dira | S B per pT bin %g-%g-%g-%g-%g-%g | all S B S/B Seff\n", ptb[0], ptb[1], ptb[2], ptb[3], ptb[4], ptb[5]);
  for (double dpi : {0.02, 0.05, 0.1, 0.15, 0.2, 0.3})
    for (double dp : {0.02, 0.03, 0.04, 0.05, 0.1})
      for (double L : {0.5, 1.0, 2.0, 3.0})
        for (double D : {0.995, 0.998, 0.999, 0.9995, 0.9997, 0.9999}) {
          std::vector<double> ns(6, 0), nb(6, 0);
          for (const auto& c : v) {
            if (c.dpi <= dpi || c.dp <= dp || c.lxy <= L || c.dira <= D || c.pt < 0.5 || c.pt >= 4.0) continue;
            int k = 0; while (c.pt >= ptb[k + 1]) ++k;
            for (int j : {k, 5}) {
              if (c.m > sLo && c.m < sHi) ns[j] += 1;
              else if ((c.m > b1 && c.m < b2) || (c.m > b3 && c.m < b4)) nb[j] += 1;
            }
          }
          printf("SCAN %.2f %.2f %.1f %.5f |", dpi, dp, L, D);
          for (int j = 0; j < 6; ++j) {
            const double B = nb[j] * wS / wB, S = ns[j] - B;
            if (j < 5) printf(" %.0f %.0f", S, B);
            else printf(" | %.0f %.0f %.3f %.1f\n", S, B, S / std::max(B, 1.), S * S / std::max(S + B, 1.));
          }
        }
}
