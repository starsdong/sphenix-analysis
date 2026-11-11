// compare_cutsets.C : storage and V0 signal of the loose vs production cut set on the same events.
// usage: root -b -q 'compare_cutsets.C+("loose_merged_dir","prod_merged_dir","tag")'
//   files: <dir>/KFP_{Ks,Lambda,D0,EventInfo}_<tag>.root
#include <TF1.h>
#include <TFile.h>
#include <TH1F.h>
#include <TSystem.h>
#include <TTree.h>

#include <cmath>
#include <cstdio>
#include <string>

namespace
{
  // offline QA cuts of the analysis note (sec. 3.8)
  double gQaDira = 0.995;
  const char *qa(const std::string &M)
  {
    static std::string s;
    s = M + Form("_DIRA>%g && ", gQaDira) + "" + M + "_decayLength_xy>0.5 && track_1_track_2_DCA<0.2 && abs(track_1_PV_DCA_xy)>0.02 && abs(track_2_PV_DCA_xy)>0.02";
    return s.c_str();
  }
  void fitSig(TTree *t, const std::string &M, double m0, double s0, double lo, double hi, double flo, double fhi, int nb, double &S, double &eS, double &mu, double &sig)
  {
    static int k = 0;
    std::string hn = Form("h%d", k++);
    t->Draw(Form("%s_mass>>%s(%d,%f,%f)", M.c_str(), hn.c_str(), nb, lo, hi), qa(M), "goff");
    TH1F *h = (TH1F *) gDirectory->Get(hn.c_str());
    TF1 f("f", "[0]*exp(-0.5*((x-[1])/[2])^2)+[3]+[4]*(x-[1])+[5]*(x-[1])^2", flo, fhi);
    f.SetParameters(h->GetMaximum(), m0, s0, h->GetBinContent(h->FindBin(flo)), 0, 0);
    f.SetParLimits(0, 0, 1e9);
    f.SetParLimits(1, m0 - 4 * s0, m0 + 4 * s0);
    f.SetParLimits(2, 0.3 * s0, 3 * s0);
    h->Fit(&f, "QRL0");
    double bw = h->GetBinWidth(1);
    mu = f.GetParameter(1);
    sig = std::fabs(f.GetParameter(2));
    S = f.GetParameter(0) * sig * std::sqrt(2 * M_PI) / bw;
    eS = f.GetParError(0) * sig * std::sqrt(2 * M_PI) / bw;
  }
}  // namespace

void compare_cutsets(const std::string &dl, const std::string &dp, const std::string &tag = "run79509_seg000-030", double qaDira = 0.995)
{
  gQaDira = qaDira;
  printf("offline QA cuts: DIRA>%g, L_xy>0.5 cm, DCA12<0.2 cm, daughter |DCA_xy|>0.02 cm\n", qaDira);
  const char *ch[3] = {"Ks", "Lambda", "D0"};
  const char *M[3] = {"K_S0", "Lambda0", "D0"};
  double nevL = 0, nevP = 0;
  {
    TFile a((dl + "/KFP_EventInfo_" + tag + ".root").c_str()), b((dp + "/KFP_EventInfo_" + tag + ".root").c_str());
    nevL = ((TTree *) a.Get("EventInfo"))->GetEntries();
    nevP = ((TTree *) b.Get("EventInfo"))->GetEntries();
  }
  printf("events: loose %.0f  prod %.0f  %s\n", nevL, nevP, nevL == nevP ? "(same sample: PASS)" : "(DIFFERENT: FAIL)");
  printf("%-7s %12s %12s %8s | %10s %10s %8s | %9s %9s\n", "channel", "cand loose", "cand prod", "ratio", "MB loose", "MB prod", "ratio", "kB/ev L", "kB/ev P");
  double totL = 0, totP = 0;
  for (int i = 0; i < 3; ++i)
  {
    std::string fl = dl + "/KFP_" + ch[i] + "_" + tag + ".root", fp = dp + "/KFP_" + ch[i] + "_" + tag + ".root";
    Long64_t szL = 0, szP = 0;
    FileStat_t st;
    if (!gSystem->GetPathInfo(fl.c_str(), st)) szL = st.fSize;
    if (!gSystem->GetPathInfo(fp.c_str(), st)) szP = st.fSize;
    TFile a(fl.c_str()), b(fp.c_str());
    TTree *ta = (TTree *) a.Get("DecayTree"), *tb = (TTree *) b.Get("DecayTree");
    double na = ta->GetEntries(), nb = tb->GetEntries();
    totL += szL;
    totP += szP;
    printf("%-7s %12.0f %12.0f %7.1f%% | %10.1f %10.1f %7.1f%% | %9.2f %9.2f\n", ch[i], na, nb, 100 * nb / na, szL / 1048576., szP / 1048576., 100. * szP / szL,
           szL / 1024. / nevL, szP / 1024. / nevP);
    if (i < 2)
    {
      double m0 = i == 0 ? 0.4976 : 1.1157, s0 = i == 0 ? 0.005 : 0.002;
      double lo = i == 0 ? 0.40 : 1.08, hi = i == 0 ? 0.60 : 1.20, flo = i == 0 ? 0.43 : 1.09, fhi = i == 0 ? 0.57 : 1.16;
      int nbin = i == 0 ? 100 : 120;
      double Sa, eSa, mua, sa, Sb, eSb, mub, sb;
      fitSig(ta, M[i], m0, s0, lo, hi, flo, fhi, nbin, Sa, eSa, mua, sa);
      fitSig(tb, M[i], m0, s0, lo, hi, flo, fhi, nbin, Sb, eSb, mub, sb);
      printf("        signal with offline QA cuts (incl. duplicates): loose %.0f+-%.0f (mu %.1f, sigma %.1f MeV)  prod %.0f+-%.0f (mu %.1f, sigma %.1f MeV)  kept %.1f%%\n",
             Sa, eSa, 1e3 * mua, 1e3 * sa, Sb, eSb, 1e3 * mub, 1e3 * sb, 100 * Sb / Sa);
    }
  }
  {
    std::string fl = dl + "/KFP_EventInfo_" + tag + ".root";
    FileStat_t st;
    gSystem->GetPathInfo(fl.c_str(), st);
    printf("EventInfo file %.1f MB (%.2f kB/ev), identical in both cut sets\n", st.fSize / 1048576., st.fSize / 1024. / nevL);
  }
  printf("TOTAL Ks+Lambda+D0: loose %.1f MB (%.2f kB/ev)  prod %.1f MB (%.2f kB/ev)  -> %.1f%% of loose (reduction x%.1f)\n", totL / 1048576.,
         totL / 1024. / nevL, totP / 1048576., totP / 1024. / nevP, 100 * totP / totL, totL / totP);
}
