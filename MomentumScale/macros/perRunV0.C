// perRunV0.C -- per-run summary of the covariance-free V0 mass histograms (condor 480608, one file per segment).
// Prints per run: events, Ks QA candidates in 0.44-0.56 (0.2<pT<4), per event, Gaussian(+pol1) mean/width of the Ks peak
// integrated and in V0 pT 0.2-0.9 / 0.9-4.0; writes merged subsets: seg 0-30 of run 79509 (sub31) and all others (sub256).
//   root -l -b -q 'macros/perRunV0.C+("data/perjob","data/perjob/merged")'
#include <TFile.h>
#include <TH2F.h>
#include <TH1D.h>
#include <TF1.h>
#include <TParameter.h>
#include <TSystem.h>
#include <TSystemDirectory.h>
#include <TList.h>
#include <TString.h>
#include <TRegexp.h>
#include <map>
#include <vector>
#include <string>
#include <cstdio>
#include <cmath>
#include <algorithm>

namespace pr {
struct Acc { double nev = 0; int nseg = 0; std::map<std::string, TH2F*> h; };
void add(Acc& a, TFile& f, const char* tag)
{
  auto* p = (TParameter<Long64_t>*)f.Get("nEventsProcessed");
  a.nev += p ? p->GetVal() : 0; a.nseg++;
  for (const char* n : {"hKs_qa", "hLa_qa", "hLaBar_qa", "hKs_loose", "hLa_loose", "hLaBar_loose"}) {
    TH2F* x = (TH2F*)f.Get(n);
    if (!x) continue;
    if (!a.h.count(n)) { a.h[n] = (TH2F*)x->Clone(Form("%s_%s", n, tag)); a.h[n]->SetDirectory(nullptr); }
    else a.h[n]->Add(x);
  }
}
void gfit(TH2F* h, double plo, double phi, double& m, double& em, double& s, double& es, double& n)
{
  TH1D* p = h->ProjectionY(Form("%s_p%.1f", h->GetName(), plo), h->GetXaxis()->FindBin(plo + 1e-6), h->GetXaxis()->FindBin(phi - 1e-6));
  p->Rebin(2);
  n = p->Integral(p->FindBin(0.44 + 1e-6), p->FindBin(0.56 - 1e-6));
  TF1 g("g", "gaus(0)+pol1(3)", 0.46, 0.54);
  g.SetParameters(p->GetMaximum(), 0.496, 0.005, 0, 0); g.SetParLimits(2, 0.001, 0.02); g.SetParLimits(1, 0.49, 0.50);
  p->Fit(&g, "QNR");
  m = 1e3 * g.GetParameter(1); em = 1e3 * g.GetParError(1); s = 1e3 * std::fabs(g.GetParameter(2)); es = 1e3 * g.GetParError(2);
}
}  // namespace pr

void perRunV0(const char* dir, const char* outDir)
{
  using namespace pr;
  std::map<int, Acc> runs;
  Acc sub31, sub256, all;
  TSystemDirectory d("d", dir);
  TList* fl = d.GetListOfFiles();
  std::vector<std::string> names;
  for (auto* o : *fl) { TString n = o->GetName(); if (n.BeginsWith("V0mass_DST") && n.EndsWith(".root")) names.push_back(n.Data()); }
  std::sort(names.begin(), names.end());
  for (auto& n : names) {
    // ...-000RRRRR-SSSSS.root
    const size_t k = n.rfind('-');
    const int seg = std::stoi(n.substr(k + 1, 5));
    const int run = std::stoi(n.substr(k - 8, 8));
    TFile f(Form("%s/%s", dir, n.c_str()));
    add(runs[run], f, Form("r%d", run));
    add(run == 79509 && seg <= 30 ? sub31 : sub256, f, run == 79509 && seg <= 30 ? "s31" : "s256");
    add(all, f, "all");
  }
  printf("files %zu\n", names.size());
  auto line = [](const char* lab, Acc& a) {
    double m, em, s, es, n, m1, em1, s1, es1, n1, m2, em2, s2, es2, n2;
    gfit(a.h["hKs_qa"], 0.2, 4.0, m, em, s, es, n);
    gfit(a.h["hKs_qa"], 0.2, 0.9, m1, em1, s1, es1, n1);
    gfit(a.h["hKs_qa"], 0.9, 4.0, m2, em2, s2, es2, n2);
    printf("RUN %s nseg %d nev %.0f nKs %.0f perEv %.4f | all %.2f %.2f %.2f %.2f | lowpt %.2f %.2f %.2f %.2f | highpt %.2f %.2f %.2f %.2f\n",
           lab, a.nseg, a.nev, n, n / a.nev, m, em, s, es, m1, em1, s1, es1, m2, em2, s2, es2);
  };
  for (auto& r : runs) line(Form("%d", r.first), r.second);
  line("sub31", sub31); line("sub256", sub256); line("all287", all);
  gSystem->mkdir(outDir, true);
  for (auto p : {std::make_pair("sub31", &sub31), std::make_pair("sub256", &sub256)}) {
    TFile o(Form("%s/V0mass_%s.root", outDir, p.first), "RECREATE");
    for (auto& h : p.second->h) h.second->Write(h.first.c_str());
    TParameter<Long64_t>("nEventsProcessed", (Long64_t)p.second->nev).Write();
    o.Close();
  }
  // per-run files (for per-run-group fits)
  for (auto& r : runs) {
    TFile o(Form("%s/V0mass_run%d.root", outDir, r.first), "RECREATE");
    for (auto& h : r.second.h) h.second->Write(h.first.c_str());
    TParameter<Long64_t>("nEventsProcessed", (Long64_t)r.second.nev).Write();
    o.Close();
  }
  printf("All done\n");
}
