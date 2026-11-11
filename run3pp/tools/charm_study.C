// charm_study.C : size / rate / mass-spectrum study of the 1% charm test (D+, Ds, Lc)
// usage: root -b -q 'charm_study.C+("<dir with files_<ch>.txt>", "<outdir for png/txt>")'
//   files_<ch>.txt : one KFP_<ch>_*.root per line (finished jobs only)
#include <TCanvas.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1F.h>
#include <TLatex.h>
#include <TParameter.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTree.h>

#include <cstdio>
#include <fstream>
#include <map>
#include <string>
#include <vector>

namespace
{
  std::string runClass(int run)
  {
    if (run < 80000) return "795xx";
    if (run < 81000) return "808xx-809xx";
    return "815xx-816xx";
  }
  struct Acc
  {
    double nev = 0, ncand = 0, bytes = 0;
    int nfiles = 0;
  };
}  // namespace

void charm_study(const std::string &dir, const std::string &out)
{
  gStyle->SetOptStat(0);
  gSystem->mkdir(out.c_str(), true);
  FILE *log = fopen((out + "/charm_study.txt").c_str(), "w");
  struct Ch
  {
    std::string key, M, label, pt3, tight;
    double lo, hi;
    std::vector<double> lines;
  };
  std::vector<Ch> chans = {
      {"Dplus", "Dplus", "D^{#pm} #rightarrow K^{#mp}#pi^{#pm}#pi^{#pm}", "track_1_pT>0.5&&track_2_pT>0.5&&track_3_pT>0.5",
       "Dplus_DIRA>0.995&&Dplus_decayLength/Dplus_decayLengthErr>3", 1.70, 2.05, {1.8697}},
      {"Ds", "D_splus", "D_{s}^{#pm} #rightarrow #phi#pi^{#pm} #rightarrow K^{+}K^{-}#pi^{#pm}", "phi_track_1_pT>0.5&&phi_track_2_pT>0.5&&track_3_pT>0.5",
       "D_splus_DIRA>0.995&&D_splus_decayLength/D_splus_decayLengthErr>3&&abs(phi_mass-1.0195)<0.008", 1.75, 2.10, {1.8697, 1.9683}},
      {"Lc", "Lambda_cplus", "#Lambda_{c}^{#pm} #rightarrow pK#pi", "track_1_pT>0.5&&track_2_pT>0.5&&track_3_pT>0.5",
       "Lambda_cplus_DIRA>0.995&&Lambda_cplus_decayLength/Lambda_cplus_decayLengthErr>2", 2.15, 2.42, {2.2865}}};

  TCanvas c("c", "c", 1500, 900);
  c.Divide(3, 2);
  int pad = 0;
  for (auto &ch : chans)
  {
    std::ifstream in(dir + "/files_" + ch.key + ".txt");
    std::string f;
    std::map<std::string, Acc> acc;
    Acc tot;
    TChain t("DecayTree");
    while (std::getline(in, f))
    {
      if (f.empty()) continue;
      // KFP_<ch>_<run>_<seg>_00000.root
      std::string b = gSystem->BaseName(f.c_str());
      // KFP_<ch>_<run>_<seg>_00000.root (per job) or KFP_<ch>_run<run>.root (merged)
      size_t pr = b.find("_run");
      int run = pr != std::string::npos ? std::stoi(b.substr(pr + 4, 8)) : std::stoi(b.substr(b.find('_', 4) + 1, 8));
      FileStat_t st;
      gSystem->GetPathInfo(f.c_str(), st);
      TFile tf(f.c_str());
      auto *np = (TParameter<Long64_t> *) tf.Get("nEventsProcessed");
      auto *tr = (TTree *) tf.Get("DecayTree");
      double nev = np ? np->GetVal() : 0, nc = tr ? tr->GetEntries() : 0;
      for (Acc *a : {&acc[runClass(run)], &tot})
      {
        a->nev += nev;
        a->ncand += nc;
        a->bytes += st.fSize;
        a->nfiles++;
      }
      if (tr) t.Add(f.c_str());
    }
    fprintf(log, "== %s: %d files, %.0f events, %.0f candidates (%.4f /event), %.1f MB, %.3f kB/event, %.0f B/candidate\n", ch.key.c_str(), tot.nfiles,
            tot.nev, tot.ncand, tot.ncand / tot.nev, tot.bytes / 1048576., tot.bytes / 1024. / tot.nev, tot.bytes / tot.ncand);
    for (auto &[k, a] : acc)
      fprintf(log, "   %-12s files %5d events %9.0f cand/event %.4f kB/event %.3f\n", k.c_str(), a.nfiles, a.nev, a.ncand / a.nev, a.bytes / 1024. / a.nev);
    fprintf(log, "   projected to 100%% (7.2e8 events, current DST area): %.1f GB\n", tot.bytes / tot.nev * 7.2e8 / 1e9);

    // mass spectra: all candidates, and offline-tightened selections
    std::string M = ch.M;
    std::vector<std::pair<std::string, std::string>> sel = {{"all (production cuts)", "1"},
                                                           {"+ daughters p_{T}>0.5", ch.pt3},
                                                           {"+ DIRA>0.995, L/#sigma_{L}", ch.pt3 + "&&" + ch.tight}};
    int col[3] = {kBlack, kBlue + 1, kRed + 1};
    std::vector<TH1F *> hs;
    for (int i = 0; i < 3; ++i)
    {
      std::string hn = Form("h_%s_%d", ch.key.c_str(), i);
      t.Draw(Form("%s_mass>>%s(70,%f,%f)", M.c_str(), hn.c_str(), ch.lo, ch.hi), sel[i].second.c_str(), "goff");
      TH1F *h = (TH1F *) gDirectory->Get(hn.c_str());
      h->SetLineColor(col[i]);
      h->SetLineWidth(2);
      hs.push_back(h);
      fprintf(log, "   selection %-28s : %8.0f candidates\n", sel[i].first.c_str(), h->GetEntries());
    }
    c.cd(++pad);
    hs[0]->SetTitle(Form("%s;mass (GeV/c^{2});candidates / 5 MeV", ch.label.c_str()));
    hs[0]->SetMinimum(0);
    hs[0]->Draw("hist");
    hs[1]->Draw("hist same");
    c.cd(pad + 3);
    hs[2]->SetTitle(Form("%s, tight offline;mass (GeV/c^{2});candidates / 5 MeV", ch.label.c_str()));
    hs[2]->SetMinimum(0);
    hs[2]->Draw("hist e");
    TLatex lt;
    lt.SetNDC();
    lt.SetTextSize(0.035);
    lt.DrawLatex(0.15, 0.85, sel[2].second.substr(0, 60).c_str());
    for (double m : ch.lines)
    {
      TLatex tl;
      tl.SetTextSize(0.03);
      tl.SetTextColor(kGray + 2);
      tl.DrawLatex(m, hs[2]->GetMaximum() * 0.05, "|");
    }
  }
  c.SaveAs((out + "/charm1pct_mass.png").c_str());
  fclose(log);
  printf("charm_study done\n");
}
