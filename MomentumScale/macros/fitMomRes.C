// fitMomRes.C -- pT-differential fit of the DDSim momentum response to the run3pp K0S / Lambda mass peaks.
//
// Model (per daughter, species mass m, true momentum p, beta = p/E):
//   p_rec = s * p * [1 + delta_sp(pT) + sigma(pT) G],   sigma(pT)/pT = sqrt((a pT)^2 + (b/beta)^2),   G ~ N(0,1)
// delta_sp = momentum bias of the full simulation (hMomShift_<sp> of the DDSim input; 0 for pions).
//
// Data   : TH2 (pT, M) of the V0 mass after the QA cuts (macro run3pp/macro/loosecov/fill_v0_mass.C), hadd-ed.
// DDSim  : trees generated with --no-mom-smear (daughter momenta at the reconstructed secondary vertex stored);
//          for each (a, b, s) the momenta are rescaled with FIXED random numbers G (smooth chi2 surface),
//          the mass and mother pT recomputed and filled into templates per mother-pT bin.
// Per pT bin the data are fitted with N * template + linear background (linear least squares, variance
// = max(data,1) + N^2 var(template)); chi2 summed over pT bins (and channels).
//
// Optional laTight = 1 (last argument): tight Lambda selection in data AND templates -- data histograms hLa_tight /
// hLaBar_tight (run3pp/macro/kappa/fill_v0_mass_kappa_v2.C: QA + DIRA > 0.9995, |DCA_xy(pi)| > 0.10, |DCA_xy(p)| > 0.03 cm),
// templates passAll && rcCosTheta > 0.9995 && |rcPiDcaXY| > 0.10 && |rcPDcaXY| > 0.03. Ks unchanged.
//
//   root -l -b -q 'macros/fitMomRes.C+("<data V0mass.root>","output/momfit/Ks_nosmear.root","output/momfit/Lambda_nosmear.root",
//                  "input/ddsim_input_run3pp.root","<outdir>","<tag>")'
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH2F.h>
#include <TF1.h>
#include <TRandom3.h>
#include <TMatrixD.h>
#include <TVectorD.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TGraphErrors.h>
#include <TStyle.h>
#include <TSystem.h>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <functional>

namespace fm {
struct LaTight { bool on = false; double dira = 0.9995, piDca = 0.10, pDca = 0.03; } gLaTight;
const double mPi = 0.13957039, mP = 0.93827209;

struct Cand { float p1[3], p2[3], g1, g2, w; };

struct Channel {
  std::string name;                    // Ks | Lambda
  double m1, m2;                       // daughter masses
  bool proton1 = false;                // daughter 1 is a proton (momentum bias)
  std::vector<double> ptE;             // mother-pT bins
  double fitLo, fitHi;                 // mass fit range
  int nMass; double mLo, mHi;          // template binning (= data binning)
  int rebin = 1;                       // rebin factor of the data mass axis
  std::vector<Cand> c;
  std::vector<TH1D*> data;             // per pT bin
};

TH1* gShiftP = nullptr;
double shiftP(double pt)
{
  if (!gShiftP) return 0;
  TAxis* ax = gShiftP->GetXaxis();
  return gShiftP->GetBinContent(std::min(std::max(ax->FindFixBin(pt), 1), ax->GetNbins()));
}

inline double resol(double a, double b, double p, double pt, double m)
{
  const double beta = p / std::sqrt(p * p + m * m);
  return std::sqrt(a * a * pt * pt + b * b / (beta * beta));
}

// fill templates for parameters (a, b, s)
void templates(const Channel& ch, double a, double b, double s, std::vector<TH1D*>& T)
{
  for (auto* h : T) h->Reset();
  const int nP = ch.ptE.size() - 1;
  for (const Cand& c : ch.c) {
    const double pt1 = std::hypot(c.p1[0], c.p1[1]), pt2 = std::hypot(c.p2[0], c.p2[1]);
    const double P1 = std::sqrt(pt1 * pt1 + c.p1[2] * c.p1[2]), P2 = std::sqrt(pt2 * pt2 + c.p2[2] * c.p2[2]);
    const double f1 = s * (1 + (ch.proton1 ? shiftP(pt1) : 0.) + resol(a, b, P1, pt1, ch.m1) * c.g1);
    const double f2 = s * (1 + resol(a, b, P2, pt2, ch.m2) * c.g2);
    const double x1 = f1 * c.p1[0], y1 = f1 * c.p1[1], z1 = f1 * c.p1[2];
    const double x2 = f2 * c.p2[0], y2 = f2 * c.p2[1], z2 = f2 * c.p2[2];
    const double e1 = std::sqrt(f1 * f1 * P1 * P1 + ch.m1 * ch.m1), e2 = std::sqrt(f2 * f2 * P2 * P2 + ch.m2 * ch.m2);
    const double px = x1 + x2, py = y1 + y2, pz = z1 + z2, E = e1 + e2;
    const double m2 = E * E - px * px - py * py - pz * pz;
    const double pt = std::hypot(px, py);
    if (pt < ch.ptE.front() || pt >= ch.ptE.back()) continue;
    const int j = int(std::upper_bound(ch.ptE.begin(), ch.ptE.end(), pt) - ch.ptE.begin()) - 1;
    if (j >= 0 && j < nP) T[j]->Fill(std::sqrt(std::max(m2, 0.)), c.w);
  }
}

// template + linear background fit in one pT bin, Poisson maximum likelihood (iteratively reweighted least squares
// with variance = model); returns the Baker-Cousins likelihood-ratio chi2 2*sum[mu - d + d ln(d/mu)], sets ndf, N.
// The DDSim template statistical error is neglected (templates have >= 50x the data statistics).
int gBkgOrder = 1;   // polynomial order of the background (1 = linear, 2 = quadratic)
double fitBin(const TH1D* d, const TH1D* t, double lo, double hi, int& ndf, double& N, double* bkg = nullptr)
{
  const int b1 = d->GetXaxis()->FindFixBin(lo + 1e-9), b2 = d->GetXaxis()->FindFixBin(hi - 1e-9);
  const double tInt = t->Integral(b1, b2), dInt = d->Integral(b1, b2);
  if (tInt <= 0 || dInt <= 0) { ndf = 0; N = 0; return 0; }
  const double m0 = 0.5 * (lo + hi);
  const int n = b2 - b1 + 1;
  // start: half of the counts in the signal template, half flat
  const int np = 2 + gBkgOrder;
  double par[4] = {0.5 * dInt / tInt, 0.5 * dInt / n, 0, 0};
  auto model = [&](int i) { const double x = d->GetBinCenter(i) - m0; return par[0] * t->GetBinContent(i) + par[1] + par[2] * x + par[3] * x * x; };
  for (int it = 0; it < 8; ++it) {
    TMatrixD A(np, np); TVectorD B(np);
    for (int i = b1; i <= b2; ++i) {
      const double mu = std::max(model(i), 0.5);
      const double x = d->GetBinCenter(i) - m0;
      const double g[4] = {t->GetBinContent(i), 1., x, x * x};
      for (int r = 0; r < np; ++r) { B(r) += g[r] * d->GetBinContent(i) / mu; for (int q = 0; q < np; ++q) A(r, q) += g[r] * g[q] / mu; }
    }
    if (std::fabs(A.Determinant()) < 1e-300) break;
    A.Invert();
    TVectorD sol = A * B;
    for (int r = 0; r < np; ++r) par[r] = sol(r);
  }
  double chi2 = 0;
  for (int i = b1; i <= b2; ++i) {
    const double mu = std::max(model(i), 1e-6), x = d->GetBinContent(i);
    chi2 += 2 * (mu - x + (x > 0 ? x * std::log(x / mu) : 0.));
  }
  ndf = n - np;
  N = par[0];
  if (bkg) { bkg[0] = par[1]; bkg[1] = par[2]; bkg[2] = m0; bkg[3] = par[3]; }
  return chi2;
}

double chi2Channel(Channel& ch, double a, double b, double s, std::vector<TH1D*>& T, int& ndf)
{
  templates(ch, a, b, s, T);
  double c2 = 0; ndf = 0;
  for (size_t j = 0; j < ch.data.size(); ++j) {
    int nd; double N;
    c2 += fitBin(ch.data[j], T[j], ch.fitLo, ch.fitHi, nd, N);
    ndf += nd;
  }
  return c2;
}

void loadSim(Channel& ch, const char* file, const std::string& d1, const std::string& d2, unsigned seed, long maxCand)
{
  TFile f(file);
  TTree* t = (TTree*)f.Get("ddsim");
  float p1[3], p2[3], w;
  bool pass;
  const char* c[3] = {"x", "y", "z"};
  for (int k = 0; k < 3; ++k) {
    t->SetBranchAddress(Form("rc%sSvP%s", d1.c_str(), c[k]), &p1[k]);
    t->SetBranchAddress(Form("rc%sSvP%s", d2.c_str(), c[k]), &p2[k]);
  }
  t->SetBranchAddress("wLevy", &w);
  t->SetBranchAddress("passAll", &pass);
  float cosT = 1, dPi = 1, dP = 1;
  const bool tight = gLaTight.on && ch.name == "Lambda";
  if (tight) { t->SetBranchAddress("rcCosTheta", &cosT); t->SetBranchAddress("rcPiDcaXY", &dPi); t->SetBranchAddress("rcPDcaXY", &dP); }
  TRandom3 r(seed);
  for (Long64_t i = 0; i < t->GetEntries() && (long)ch.c.size() < maxCand; ++i) {
    t->GetEntry(i);
    if (!pass) continue;
    if (tight && !(cosT > gLaTight.dira && std::fabs(dPi) > gLaTight.piDca && std::fabs(dP) > gLaTight.pDca)) continue;
    Cand x;
    for (int k = 0; k < 3; ++k) { x.p1[k] = p1[k]; x.p2[k] = p2[k]; }
    x.g1 = r.Gaus(); x.g2 = r.Gaus(); x.w = w;
    ch.c.push_back(x);
  }
  std::cout << ch.name << ": " << ch.c.size() << " DDSim candidates (passAll" << (tight ? " + tight Lambda cuts" : "") << ") from " << file << std::endl;
}

void loadData(Channel& ch, TFile& f, const std::vector<std::string>& hn)
{
  TH2* h = nullptr;
  for (auto& n : hn) {
    TH2* x = (TH2*)f.Get(n.c_str());
    if (!x) { std::cerr << "missing " << n << std::endl; continue; }
    if (!h) { h = (TH2*)x->Clone(("sum_" + ch.name).c_str()); h->SetDirectory(nullptr); }
    else h->Add(x);
  }
  if (ch.rebin > 1) h->RebinY(ch.rebin);
  ch.nMass = h->GetNbinsY(); ch.mLo = h->GetYaxis()->GetXmin(); ch.mHi = h->GetYaxis()->GetXmax();
  for (size_t j = 0; j + 1 < ch.ptE.size(); ++j) {
    const int b1 = h->GetXaxis()->FindBin(ch.ptE[j] + 1e-6), b2 = h->GetXaxis()->FindBin(ch.ptE[j + 1] - 1e-6);
    TH1D* p = h->ProjectionY(Form("data_%s_pt%zu", ch.name.c_str(), j), b1, b2);
    p->SetDirectory(nullptr);
    p->SetTitle("");
    ch.data.push_back(p);
  }
}
}  // namespace fm

void fitMomRes(const char* dataFile, const char* simKs, const char* simLa, const char* inputFile, const char* outDir,
               const char* tag = "cov", long maxCand = 2000000, int bkgOrder = 1, double laLo = 1.095, double laHi = 1.14,
               double sLo = 0.984, double sHi = 1.000, int laTight = 0)
{
  using namespace fm;
  gSystem->mkdir(outDir, true);
  gStyle->SetOptStat(0);
  { TFile fi(inputFile); TH1* h = (TH1*)fi.Get("hMomShift_p"); if (h) { gShiftP = (TH1*)h->Clone("shiftP"); gShiftP->SetDirectory(nullptr); } }

  Channel ks{"Ks", mPi, mPi, false, {0.2, 0.5, 0.7, 0.9, 1.2, 1.6, 2.2, 4.0}, 0.44, 0.56};
  gBkgOrder = bkgOrder;
  Channel la{"Lambda", mP, mPi, true, {0.5, 0.9, 1.2, 1.6, 2.4, 4.0}, laLo, laHi};
  ks.rebin = 2;   // 2 MeV/c^2 bins
  la.rebin = 2;   // 0.5 MeV/c^2 bins
  TFile fd(dataFile);
  loadData(ks, fd, {"hKs_qa"});
  gLaTight.on = laTight != 0;
  if (gLaTight.on) std::cout << "TIGHT Lambda: DIRA > " << gLaTight.dira << ", |DCAxy(pi)| > " << gLaTight.piDca << ", |DCAxy(p)| > " << gLaTight.pDca << std::endl;
  loadData(la, fd, gLaTight.on ? std::vector<std::string>{"hLa_tight", "hLaBar_tight"} : std::vector<std::string>{"hLa_qa", "hLaBar_qa"});
  loadSim(ks, simKs, "PiP", "PiM", 1, maxCand);
  loadSim(la, simLa, "P", "Pi", 2, maxCand);
  for (Channel* ch : {&ks, &la}) {
    std::cout << ch->name << " data entries per pT bin (fit range):";
    for (auto* h : ch->data) std::cout << " " << h->Integral(h->FindBin(ch->fitLo), h->FindBin(ch->fitHi));
    std::cout << std::endl;
  }
  auto mkT = [](Channel& ch) {
    std::vector<TH1D*> T;
    for (size_t j = 0; j + 1 < ch.ptE.size(); ++j) {
      TH1D* h = new TH1D(Form("tmpl_%s_%zu", ch.name.c_str(), j), "", ch.nMass, ch.mLo, ch.mHi);
      h->Sumw2(); h->SetDirectory(nullptr);
      T.push_back(h);
    }
    return T;
  };
  std::vector<TH1D*> Tk = mkT(ks), Tl = mkT(la);

  // ---------------- grid scan ----------------
  std::vector<double> A, B, S;
  for (double a = 0.0; a <= 0.08001; a += 0.0025) A.push_back(a);
  for (double b = 0.004; b <= 0.04001; b += 0.002) B.push_back(b);
  for (double s = sLo; s <= sHi + 1e-5; s += 0.001) S.push_back(s);   // default 0.984-1.000 (17 points)
  const int nA = A.size(), nB = B.size(), nS = S.size();
  std::vector<double> cK(nA * nB * nS), cL(nA * nB * nS);
  int ndfK = 0, ndfL = 0;
  for (int i = 0; i < nA; ++i) {
    for (int j = 0; j < nB; ++j)
      for (int k = 0; k < nS; ++k) {
        const int id = (i * nB + j) * nS + k;
        cK[id] = chi2Channel(ks, A[i], B[j], S[k], Tk, ndfK);
        cL[id] = chi2Channel(la, A[i], B[j], S[k], Tl, ndfL);
      }
    std::cout << "  grid a = " << A[i] << " done" << std::endl;
  }
  // ---------------- results ----------------
  std::ofstream log(Form("%s/fitMomRes_%s.txt", outDir, tag));
  TFile out(Form("%s/fitMomRes_%s.root", outDir, tag), "RECREATE");
  struct Best { double a, b, s, chi2; int ndf; };
  auto analyse = [&](const char* name, std::function<double(int)> c2, int ndf) {
    int best = 0;
    for (int id = 0; id < nA * nB * nS; ++id) if (c2(id) < c2(best)) best = id;
    const int bi = best / (nB * nS), bj = (best / nS) % nB, bk = best % nS;
    // profiles (minimum over the other parameters)
    std::vector<double> pa(nA, 1e30), pb(nB, 1e30), ps(nS, 1e30);
    TH2D* map = new TH2D(Form("chi2_ab_%s", name), Form("%s: min_{s} #chi^{2} - #chi^{2}_{min};a (1/(GeV/c));b", name),
                         nA, A.front() - 0.00125, A.back() + 0.00125, nB, B.front() - 0.001, B.back() + 0.001);
    for (int i = 0; i < nA; ++i)
      for (int j = 0; j < nB; ++j)
        for (int k = 0; k < nS; ++k) {
          const double v = c2((i * nB + j) * nS + k);
          pa[i] = std::min(pa[i], v); pb[j] = std::min(pb[j], v); ps[k] = std::min(ps[k], v);
          const double cur = map->GetBinContent(i + 1, j + 1);
          if (cur == 0 || v - c2(best) < cur) map->SetBinContent(i + 1, j + 1, std::max(v - c2(best), 1e-6));
        }
    map->Write();
    // 1-sigma interval from the profile (Delta chi2 = 1), linear interpolation on the grid
    auto interval = [](const std::vector<double>& x, const std::vector<double>& p, double& lo, double& hi) {
      const double mn = *std::min_element(p.begin(), p.end());
      lo = x.front(); hi = x.back();
      int im = std::min_element(p.begin(), p.end()) - p.begin();
      for (int i = im; i > 0; --i) if (p[i - 1] - mn > 1) { lo = x[i - 1] + (x[i] - x[i - 1]) * (1 + mn - p[i - 1]) / (p[i] - p[i - 1]); break; }
      for (int i = im; i + 1 < (int)x.size(); ++i) if (p[i + 1] - mn > 1) { hi = x[i] + (x[i + 1] - x[i]) * (1 + mn - p[i]) / (p[i + 1] - p[i]); break; }
    };
    double aLo, aHi, bLo, bHi, sLo, sHi;
    interval(A, pa, aLo, aHi); interval(B, pb, bLo, bHi); interval(S, ps, sLo, sHi);
    char buf[600];
    snprintf(buf, sizeof buf, "FIT %s a %.4f %.4f %.4f b %.4f %.4f %.4f s %.4f %.4f %.4f chi2 %.1f ndf %d edge_a %d edge_b %d edge_s %d",
             name, A[bi], aLo, aHi, B[bj], bLo, bHi, S[bk], sLo, sHi, c2(best), ndf - 3, int(bi == 0 || bi == nA - 1),
             int(bj == 0 || bj == nB - 1), int(bk == 0 || bk == nS - 1));
    std::cout << buf << std::endl; log << buf << "\n";
    return Best{A[bi], B[bj], S[bk], c2(best), ndf - 3};
  };
  Best bK = analyse("Ks", [&](int id) { return cK[id]; }, ndfK);
  Best bL = analyse("Lambda", [&](int id) { return cL[id]; }, ndfL);
  Best bJ = analyse("joint", [&](int id) { return cK[id] + cL[id]; }, ndfK + ndfL);
  // cross-check chi2 of each channel at the other's / joint best point
  int nd;
  auto at = [&](Channel& ch, std::vector<TH1D*>& T, const Best& p) { return chi2Channel(ch, p.a, p.b, p.s, T, nd); };
  char buf[400];
  snprintf(buf, sizeof buf, "CROSS chi2(Lambda | Ks best) %.1f / %d ; chi2(Lambda | own best) %.1f ; chi2(Ks | joint) %.1f ; chi2(Lambda | joint) %.1f",
           at(la, Tl, bK), ndfL, bL.chi2, at(ks, Tk, bJ), at(la, Tl, bJ));
  std::cout << buf << std::endl; log << buf << "\n";

  // ---------------- per-bin Gaussian widths (display) + mass plots at the Ks best point ----------------
  for (Channel* ch : {&ks, &la}) {
    std::vector<TH1D*>& T = ch == &ks ? Tk : Tl;
    templates(*ch, bK.a, bK.b, bK.s, T);
    const int nP = ch->data.size();
    TCanvas c(Form("c_%s", ch->name.c_str()), "", 1500, 300 * ((nP + 3) / 4));
    c.Divide(4, (nP + 3) / 4);
    TGraphErrors* gD = new TGraphErrors(); gD->SetName(Form("gWidthData_%s", ch->name.c_str()));
    TGraphErrors* gS = new TGraphErrors(); gS->SetName(Form("gWidthSim_%s", ch->name.c_str()));
    for (int j = 0; j < nP; ++j) {
      c.cd(j + 1);
      gPad->SetLeftMargin(0.15); gPad->SetTopMargin(0.12);
      int ndf; double N, bkg[4];
      const double chi2 = fitBin(ch->data[j], T[j], ch->fitLo, ch->fitHi, ndf, N, bkg);
      TH1D* d = ch->data[j];
      d->GetXaxis()->SetRangeUser(ch->fitLo, ch->fitHi);
      d->SetMarkerStyle(20); d->SetMarkerSize(0.6);
      d->GetXaxis()->SetTitle("M (GeV/c^{2})");
      d->SetMinimum(0);
      d->Draw("E");
      TH1D* m = (TH1D*)T[j]->Clone(Form("model_%s_%d", ch->name.c_str(), j));
      for (int i = 1; i <= m->GetNbinsX(); ++i) m->SetBinContent(i, N * T[j]->GetBinContent(i) + bkg[0] + bkg[1] * (m->GetBinCenter(i) - bkg[2]) + bkg[3] * std::pow(m->GetBinCenter(i) - bkg[2], 2));
      m->SetLineColor(kRed + 1); m->SetLineWidth(2);
      m->GetXaxis()->SetRangeUser(ch->fitLo, ch->fitHi);
      m->Draw("hist same");
      TLatex l; l.SetNDC(); l.SetTextSize(0.06);
      l.DrawLatex(0.17, 0.92, Form("%s  %.1f<p_{T}<%.1f GeV/c", ch->name == "Ks" ? "K^{0}_{S}" : "#Lambda+#bar{#Lambda}", ch->ptE[j], ch->ptE[j + 1]));
      l.SetTextSize(0.055);
      l.DrawLatex(0.60, 0.82, Form("#chi^{2}/ndf = %.0f/%d", chi2, ndf));
      // Gaussian widths (data: gaus + pol1; sim template: gaus) for display
      const double mPDG = ch->name == "Ks" ? 0.4976 : 1.1157;
      TF1 fdd("fdd", "gaus(0)+pol1(3)", ch->fitLo, ch->fitHi);
      fdd.SetParameters(d->GetMaximum(), mPDG, ch->name == "Ks" ? 0.005 : 0.002, 0, 0);
      fdd.SetParLimits(2, 0.0003, 0.03);
      d->Fit(&fdd, "QNR");
      TF1 fss("fss", "gaus", ch->fitLo, ch->fitHi);
      fss.SetParameters(T[j]->GetMaximum(), T[j]->GetBinCenter(T[j]->GetMaximumBin()), ch->name == "Ks" ? 0.005 : 0.002);
      T[j]->Fit(&fss, "QNR");
      const double x = 0.5 * (ch->ptE[j] + ch->ptE[j + 1]);
      gD->SetPoint(j, x, 1e3 * std::fabs(fdd.GetParameter(2))); gD->SetPointError(j, 0.5 * (ch->ptE[j + 1] - ch->ptE[j]), 1e3 * fdd.GetParError(2));
      gS->SetPoint(j, x, 1e3 * std::fabs(fss.GetParameter(2))); gS->SetPointError(j, 0, 1e3 * fss.GetParError(2));
      snprintf(buf, sizeof buf, "BIN %s %.2f %.2f chi2 %.1f ndf %d N %.3g gaus_data %.3f %.3f gaus_sim %.3f", ch->name.c_str(), ch->ptE[j],
               ch->ptE[j + 1], chi2, ndf, N, 1e3 * std::fabs(fdd.GetParameter(2)), 1e3 * fdd.GetParError(2), 1e3 * std::fabs(fss.GetParameter(2)));
      std::cout << buf << std::endl; log << buf << "\n";
      d->Write(); m->Write();
    }
    c.SaveAs(Form("%s/fitMomRes_mass_%s_%s.png", outDir, ch->name.c_str(), tag));
    c.SaveAs(Form("%s/fitMomRes_mass_%s_%s.pdf", outDir, ch->name.c_str(), tag));
    gD->Write(); gS->Write();
  }
  out.Close();
  std::cout << "All done" << std::endl;
}
