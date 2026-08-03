#include "/Users/wuwenya/Documents/code/femto_flow/CATS/CATS_install/install/include/CATS.h"
#include "/Users/wuwenya/Documents/code/femto_flow/CATS/CATS_install/install/include/CATStools.h"
#include "/Users/wuwenya/Documents/code/femto_flow/CATS/CATS_install/install/include/CATSconstants.h"
#include "/Users/wuwenya/Documents/code/femto_flow/CATS/CATS_install/install/include/DLM_Ck.h"
#include "/Users/wuwenya/Documents/code/femto_flow/CATS/CATS_install/install/include/DLM_CkDecomposition.h"
#include "/Users/wuwenya/Documents/code/femto_flow/CATS/CATS_install/install/include/DLM_Histo.h"
#include "/Users/wuwenya/Documents/code/femto_flow/CATS/CATS_install/install/include/DLM_Random.h"
#include "/Users/wuwenya/Documents/code/femto_flow/CATS/CATS_install/install/include/DLM_Source.h"

#include "TCanvas.h"
#include "TFile.h"
#include "TLegend.h"
#include "TH1.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TLine.h"
#include "TNtuple.h"
#include "TPad.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TGraph.h"

#include <algorithm>
#include <cmath>
#include <iostream>

namespace {

const TString kLambdaFile = "../Files/lam/kStarDepLambda_cent0to10_pion_PbPb.root";
const TString kLambdaHist = "hPurityKstar_mT0.40to0.50_cent0to10";
const TString kSmearFile = "../Files/resolution_rebin4.root";
const TString kSmearHist = "resolution";
const TString kPRFile = "../Files/ForMax_pi_piReso_noOmega.root";
const TString kRRFile = "../Files/ForMax_piReso_piReso_noOmega.root";

TH2F* LoadSmear()
{
  TFile f(kSmearFile, "READ");
  TH2F* h = dynamic_cast<TH2F*>(f.Get(kSmearHist));
  if(!h) {
    std::cout << "ERROR: cannot load smearing matrix\n";
    return nullptr;
  }
  TH2F* copy = dynamic_cast<TH2F*>(h->Clone("hMomentumResolution_copy_fixedR7"));
  copy->SetDirectory(nullptr);
  return copy;
}

DLM_Histo<double>* LoadFeedDown()
{
  TFile f(kLambdaFile, "READ");
  TH1* hLambda = dynamic_cast<TH1*>(f.Get(kLambdaHist));
  if(!hLambda) {
    std::cout << "ERROR: cannot load lambda hist\n";
    return nullptr;
  }

  const int nBins = hLambda->GetNbinsX();
  DLM_Histo<double>* hFeed = new DLM_Histo<double>();
  hFeed->SetUp(1);
  hFeed->SetUp(0, nBins,
               hLambda->GetXaxis()->GetBinLowEdge(1)*1000.0,
               hLambda->GetXaxis()->GetBinUpEdge(nBins)*1000.0);
  hFeed->Initialize();

  const double weakScale = (1.0-0.12)*(1.0-0.12);
  unsigned bin[1];
  for(int i=0; i<nBins; ++i) {
    const double lambdaMain = hLambda->GetBinContent(i+1)*weakScale;
    bin[0] = i;
    hFeed->SetBinContent(bin, 1.0-lambdaMain);
    hFeed->SetBinError(bin, hLambda->GetBinError(i+1)*weakScale);
  }
  return hFeed;
}

void SetPipi(CATS& kitty)
{
  kitty.SetMomBins(500, 0.0, 400.0);
  kitty.SetNumChannels(1);
  kitty.SetNumPW(0, 0);
  kitty.SetSpin(0, 0);
  kitty.SetChannelWeight(0, 1.0);
  kitty.SetQ1Q2(1);
  kitty.SetPdgId(211, 211);
  kitty.SetRedMass(0.5*Mass_Pic);
  kitty.SetNotifications(CATS::nWarning);
  kitty.SetNotifications(CATS::nError);
}

void SetGauss(CATS& kitty, const double radius)
{
  CATSparameters pars(CATSparameters::tSource, 1, true);
  pars.SetParameter(0, radius);
  kitty.SetAnaSource(GaussSource, pars);
  kitty.SetUseAnalyticSource(true);
  kitty.SetAutoNormSource(false);
}

DLM_CleverMcLevyResoTM* SetReso(CATS& kitty,
                                const double radius,
                                const double alpha,
                                const double fracReso,
                                const double tauReso,
                                const double cutoff)
{
  DLM_CleverMcLevyResoTM* source = new DLM_CleverMcLevyResoTM();
  source->InitStability(1, alpha-1e-6, alpha+1e-6);
  source->InitScale(160, 4.0, 30.0);
  source->InitRad(400, 1e-4, 128.0);
  source->InitType(2);
  source->SetUpReso(0, fracReso);
  source->SetUpReso(1, fracReso);

  Float_t k_D=0, fP1=0, fP2=0, fM1=0, fM2=0, Tau1=0, Tau2=0;
  Float_t AngleRcP1=0, AngleRcP2=0, AngleP1P2=0;
  DLM_Random ranGen(11);

  TFile filePR(kPRFile, "READ");
  TNtuple* tuplePR = dynamic_cast<TNtuple*>(filePR.Get("InfoTuple_ClosePairs"));
  if(tuplePR) {
    tuplePR->SetBranchAddress("k_D", &k_D);
    tuplePR->SetBranchAddress("P1", &fP1);
    tuplePR->SetBranchAddress("P2", &fP2);
    tuplePR->SetBranchAddress("M1", &fM1);
    tuplePR->SetBranchAddress("M2", &fM2);
    tuplePR->SetBranchAddress("Tau1", &Tau1);
    tuplePR->SetBranchAddress("Tau2", &Tau2);
    tuplePR->SetBranchAddress("AngleRcP1", &AngleRcP1);
    tuplePR->SetBranchAddress("AngleRcP2", &AngleRcP2);
    tuplePR->SetBranchAddress("AngleP1P2", &AngleP1P2);
    unsigned added = 0;
    for(Long64_t i=0; i<tuplePR->GetEntries(); ++i) {
      tuplePR->GetEntry(i);
      if(k_D > cutoff) continue;
      if(fP2 <= 0.f || fM2 <= 0.f || tauReso <= 0.0) continue;
      const double shift = ranGen.Exponential(fM2/(fP2*tauReso));
      source->AddBGT_PR(shift, -std::cos(AngleRcP2));
      source->AddBGT_RP(shift,  std::cos(AngleRcP2));
      ++added;
    }
    std::cout << "FixedR loaded PR/RP kernels: " << added << "\n";
  }

  TFile fileRR(kRRFile, "READ");
  TNtuple* tupleRR = dynamic_cast<TNtuple*>(fileRR.Get("InfoTuple_ClosePairs"));
  if(tupleRR) {
    tupleRR->SetBranchAddress("k_D", &k_D);
    tupleRR->SetBranchAddress("P1", &fP1);
    tupleRR->SetBranchAddress("P2", &fP2);
    tupleRR->SetBranchAddress("M1", &fM1);
    tupleRR->SetBranchAddress("M2", &fM2);
    tupleRR->SetBranchAddress("Tau1", &Tau1);
    tupleRR->SetBranchAddress("Tau2", &Tau2);
    tupleRR->SetBranchAddress("AngleRcP1", &AngleRcP1);
    tupleRR->SetBranchAddress("AngleRcP2", &AngleRcP2);
    tupleRR->SetBranchAddress("AngleP1P2", &AngleP1P2);
    unsigned added = 0;
    for(Long64_t i=0; i<tupleRR->GetEntries(); ++i) {
      tupleRR->GetEntry(i);
      if(k_D > cutoff) continue;
      if(fM1 <= 0.f || fM2 <= 0.f) continue;
      if(fP1 <= 0.f || fP2 <= 0.f || tauReso <= 0.0) continue;
      const double shift1 = ranGen.Exponential(fM1/(fP1*tauReso));
      const double shift2 = ranGen.Exponential(fM2/(fP2*tauReso));
      source->AddBGT_RR(shift1, std::cos(AngleRcP1),
                        shift2, std::cos(AngleRcP2),
                        std::cos(AngleP1P2));
      ++added;
    }
    std::cout << "FixedR loaded RR kernels: " << added << "\n";
  }

  source->InitNumMcIter(1000000);
  kitty.SetAnaSource(CatsSourceForwarder, source, 2);
  kitty.SetAnaSource(0, radius);
  kitty.SetAnaSource(1, alpha);
  kitty.SetUseAnalyticSource(true);
  kitty.SetAutoNormSource(false);
  return source;
}

TH1D* BuildFixedRModel(const char* name,
                       const bool useResonance,
                       const double radius,
                       const TH1* hTemplate)
{
  CATS kitty;
  SetPipi(kitty);

  DLM_CleverMcLevyResoTM* reso = nullptr;
  if(useResonance) reso = SetReso(kitty, radius, 2.0, 0.72, 1.5, 400.0);
  else SetGauss(kitty, radius);
  kitty.KillTheCat();

  DLM_Ck ck(kitty.GetNumSourcePars(), 0, kitty);
  ck.Update(true);

  TH2F* hSmear = LoadSmear();
  DLM_Histo<double>* hFeed = LoadFeedDown();
  DLM_CkDecomposition dec(name, 1, ck, hSmear);
  if(hFeed) dec.AddContribution(0, *hFeed, DLM_CkDecomposition::cFeedDown);
  dec.Update(true, true);

  TH1D* h = dynamic_cast<TH1D*>(hTemplate->Clone(name));
  if(!h) h = new TH1D(name, name, 100, 0.0, 200.0);
  h->SetDirectory(nullptr);
  h->Reset("ICES");
  for(int iBin=1; iBin<=h->GetNbinsX(); ++iBin) {
    const double k = h->GetXaxis()->GetBinCenter(iBin);
    h->SetBinContent(iBin, dec.EvalCk(k));
    h->SetBinError(iBin, 0.0);
  }

  delete hFeed;
  delete hSmear;
  delete reso;
  return h;
}

TH1D* MakeRatio(const TH1* hNum, const TH1* hDen, const char* name)
{
  TH1D* hRatio = dynamic_cast<TH1D*>(hNum->Clone(name));
  hRatio->SetDirectory(nullptr);
  hRatio->Reset("ICES");
  for(int iBin=1; iBin<=hRatio->GetNbinsX(); ++iBin) {
    const double den = hDen->GetBinContent(iBin);
    if(den==0.) continue;
    hRatio->SetBinContent(iBin, hNum->GetBinContent(iBin)/den);
  }
  return hRatio;
}

TGraph* MakeGraphInRange(const TH1* h, const char* name, const double xmin, const double xmax)
{
  TGraph* g = new TGraph();
  g->SetName(name);
  int iPoint = 0;
  for(int iBin=1; iBin<=h->GetNbinsX(); ++iBin) {
    const double x = h->GetXaxis()->GetBinCenter(iBin);
    if(x < xmin || x > xmax) continue;
    const double y = h->GetBinContent(iBin);
    if(y <= 0.0) continue;
    g->SetPoint(iPoint++, x, y);
  }
  return g;
}

}

void compare_fixedR7_gauss_reso_lowK_noData()
{
  gStyle->SetOptStat(0);

  const double fixedR = 7.0;
  const TString outDir = "/Users/wuwenya/Documents/code/femto_flow/CATS/CATS_install/CATS_Pion_PbPb/CMake/out/Aug";
  gSystem->mkdir(outDir, true);

  TH1D* hTemplate = new TH1D("hTheoryTemplate_lowK", "", 1200, 0.0, 120.0);
  hTemplate->SetDirectory(nullptr);

  TH1D* hGauss = BuildFixedRModel("hFixedR7_pureGauss", false, fixedR, hTemplate);
  TH1D* hReso = BuildFixedRModel("hFixedR7_gaussReso_alpha2", true, fixedR, hTemplate);
  TH1D* hRatio = MakeRatio(hReso, hGauss, "hRatio_reso_over_gauss_fixedR7");

  double maxAbs = 0.0;
  double maxRel = 0.0;
  double maxK = 0.0;
  for(int iBin=1; iBin<=hGauss->GetNbinsX(); ++iBin) {
    const double k = hGauss->GetXaxis()->GetBinCenter(iBin);
    if(k < 2.0 || k > 200.0) continue;
    const double g = hGauss->GetBinContent(iBin);
    const double r = hReso->GetBinContent(iBin);
    const double diff = r - g;
    const double rel = g!=0.0 ? diff/g : 0.0;
    if(std::fabs(diff) > std::fabs(maxAbs)) {
      maxAbs = diff;
      maxRel = rel;
      maxK = k;
    }
  }

  std::cout << "Fixed R=7 fm, alpha=2 comparison\n";
  std::cout << "  maxAbs(reso-gauss) = " << maxAbs
            << " at k*=" << maxK << " MeV/c"
            << " maxRel=" << maxRel << "\n";
  for(double k : {6., 10., 14., 18., 22., 30., 40., 60., 80., 100.}) {
    const int bin = hGauss->FindBin(k);
    const double g = hGauss->GetBinContent(bin);
    const double r = hReso->GetBinContent(bin);
    std::cout << Form("  k=%6.1f  Gauss=%.8f  Reso=%.8f  diff=%+.8f  ratio=%.8f\n",
                      k, g, r, r-g, g!=0.0 ? r/g : 0.0);
  }

  hGauss->SetLineColor(kMagenta+1);
  hGauss->SetLineWidth(3);
  hReso->SetLineColor(kGreen+2);
  hReso->SetLineWidth(3);
  hReso->SetLineStyle(2);
  hRatio->SetLineColor(kGreen+2);
  hRatio->SetLineWidth(2);
  TGraph* gGauss = MakeGraphInRange(hGauss, "gFixedR7_pureGauss", 0.0, 120.0);
  TGraph* gReso = MakeGraphInRange(hReso, "gFixedR7_gaussReso_alpha2", 0.0, 120.0);
  gGauss->SetLineColor(kMagenta+1);
  gGauss->SetLineWidth(3);
  gReso->SetLineColor(kGreen+2);
  gReso->SetLineWidth(3);
  gReso->SetLineStyle(2);

  TCanvas c("c_fixedR7_gauss_reso_lowK", "fixed R=7 Gauss vs resonance low k*", 900, 800);
  TPad p1("p1", "p1", 0.0, 0.28, 1.0, 1.0);
  TPad p2("p2", "p2", 0.0, 0.0, 1.0, 0.28);
  p1.SetBottomMargin(0.02);
  p2.SetTopMargin(0.04);
  p2.SetBottomMargin(0.32);
  p1.Draw();
  p2.Draw();

  p1.cd();
  TH1D* hFrame = dynamic_cast<TH1D*>(hTemplate->Clone("hFrame_fixedR7_lowK"));
  hFrame->SetDirectory(nullptr);
  hFrame->Reset("ICES");
  hFrame->GetXaxis()->SetRangeUser(0.0, 60.0);
  hFrame->GetYaxis()->SetRangeUser(0.94, 1.36);
  hFrame->GetYaxis()->SetTitle("C(#it{k*})");
  hFrame->GetXaxis()->SetLabelSize(0.0);
  hFrame->SetTitle("Fixed R=7 fm: pure Gaussian vs Gaussian+resonance, low k*");
  hFrame->Draw("axis");
  gGauss->Draw("L same");
  gReso->Draw("L same");

  TLegend leg(0.48, 0.62, 0.88, 0.88);
  leg.SetBorderSize(0);
  leg.SetFillStyle(0);
  leg.AddEntry(gGauss, "pure Gauss, R=7 fm", "l");
  leg.AddEntry(gReso, "Gauss+ForMax Resonance, R=7 fm", "l");
  leg.Draw();

  p2.cd();
  hRatio->GetXaxis()->SetRangeUser(0.0, 60.0);
  hRatio->GetYaxis()->SetRangeUser(0.985, 1.005);
  hRatio->GetYaxis()->SetTitle("Reso/Gauss");
  hRatio->GetXaxis()->SetTitle("#it{k*} (MeV/#it{c})");
  hRatio->GetYaxis()->SetNdivisions(505);
  hRatio->GetYaxis()->SetTitleSize(0.10);
  hRatio->GetYaxis()->SetTitleOffset(0.45);
  hRatio->GetYaxis()->SetLabelSize(0.08);
  hRatio->GetXaxis()->SetTitleSize(0.10);
  hRatio->GetXaxis()->SetLabelSize(0.08);
  hRatio->Draw("hist");
  TLine one(0.0, 1.0, 60.0, 1.0);
  one.SetLineStyle(2);
  one.SetLineColor(kGray+2);
  one.Draw("same");

  const TString outPng = outDir + "/compare_fixedR7_gauss_vs_ForMaxReso_alpha2_lowK.png";
  const TString outPdf = outDir + "/compare_fixedR7_gauss_vs_ForMaxReso_alpha2_lowK.pdf";
  const TString outRoot = outDir + "/compare_fixedR7_gauss_vs_ForMaxReso_alpha2_lowK.root";
  c.Print(outPng);
  c.Print(outPdf);

  TFile out(outRoot, "RECREATE");
  hTemplate->Write();
  hGauss->Write();
  hReso->Write();
  hRatio->Write();
  hFrame->Write();
  gGauss->Write();
  gReso->Write();
  c.Write();
  out.Close();

  std::cout << "Saved fixed-R comparison:\n"
            << "  " << outPng << "\n"
            << "  " << outRoot << "\n";
}
