#include "Basics.h"
#include "ExtendedCk.h"

#include "/Users/wuwenya/Documents/code/femto_flow/CATS/CATS_install/install/include/CATS.h"
#include "/Users/wuwenya/Documents/code/femto_flow/CATS/CATS_install/install/include/CATStools.h"
#include "/Users/wuwenya/Documents/code/femto_flow/CATS/CATS_install/install/include/CATSconstants.h"
#include "/Users/wuwenya/Documents/code/femto_flow/CATS/CATS_install/install/include/DLM_Ck.h"
#include "/Users/wuwenya/Documents/code/femto_flow/CATS/CATS_install/install/include/DLM_CkDecomposition.h"
#include "/Users/wuwenya/Documents/code/femto_flow/CATS/CATS_install/install/include/DLM_Histo.h"
#include "/Users/wuwenya/Documents/code/femto_flow/CATS/CATS_install/install/include/DLM_Source.h"

#include "TFile.h"
#include "TF1.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "TGraph.h"
#include "TGraphAsymmErrors.h"
#include "TGraphErrors.h"
#include "TH1.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TMath.h"
#include "TNtuple.h"
#include "TROOT.h"
#include "TString.h"

#include <cmath>
#include <iostream>
#include <vector>

namespace PipiIO {
  TString DataFileName = "../Files/cf/cf_data_cent10to30_allTOF_pippip_cpr0.028.root";
  TString DataHistoName = "cf_same_proj_mt4_cent1_scaled";

  TString LambdaFileName = "../Files/lam/kStarDepLambda_cent10to30_pion_PbPb.root";
  TString LambdaHistoName = "hPurityKstar_mT0.60to0.70_cent10to30";

  const TString SmearingFileName = "../Files/resolution_rebin4.root";
  const TString SmearingHistoName = "resolution";

  // ForMax resonance transfer model input. The PbPb PRF files are kept here
  // for reference:
  //   ../Files/EPOSPbPb_pi_piReso_noOmega_PRF.root
  //   ../Files/EPOSPbPb_piReso_piReso_noOmega_PRF.root
  const TString LegacyPRFileName = "../Files/ForMax_pi_piReso_noOmega.root";
  const TString LegacyRRFileName = "../Files/ForMax_piReso_piReso_noOmega.root";

  TString OutputFitGaussOnly1D = "./out/Aug/fit_results_gaussOnly_ForMax_cent10to30_mT4.root";
  TString OutputFitCoreReso1D = "./out/Aug/fit_results_CleverMcLevyResoTM_ForMax_fixedAlpha2p0_cent10to30_mT4.root";

  const unsigned CatsNumMomBins = 500;
  const double CatsKMin = 0.0;
  const double CatsKMax = 400.0;
  const int FitNpx = 500;
  const double FitMin = 2.0;
  const double FitMax = 200.0;
  const double BaselineFitMin = 200.0;
  const double BaselineFitMax = 400.0;
}

static CATS* gKittyPiPi = nullptr;
static DLM_Ck* gCkPiPi = nullptr;
static DLM_CleverMcLevyResoTM* gCleverMcLevyResoTM_pipi = nullptr;
static DLM_CkDecomposition* gFemto_pipi_fitter = nullptr;

double fracWeakDecay = 0.12;
const double fracShortLivedResonance = 0.72;
double sourceAlpha = 2.0;

namespace {

void CleanupPipiTheory()
{
  if(gFemto_pipi_fitter){ delete gFemto_pipi_fitter; gFemto_pipi_fitter = nullptr; }
  if(gCkPiPi){ delete gCkPiPi; gCkPiPi = nullptr; }
  if(gKittyPiPi){ delete gKittyPiPi; gKittyPiPi = nullptr; }
  if(gCleverMcLevyResoTM_pipi){
    delete gCleverMcLevyResoTM_pipi;
    gCleverMcLevyResoTM_pipi = nullptr;
  }
}

TH2F* GetSmearMatrix(TString filename, TString histoname)
{
  TFile* file = TFile::Open(filename, "READ");
  if(!file || file->IsZombie()){
    std::cout << "ERROR GetSmearMatrix: cannot open " << filename << "\n";
    if(file){ file->Close(); delete file; }
    return nullptr;
  }

  TH2F* histo = dynamic_cast<TH2F*>(file->Get(histoname));
  if(!histo){
    std::cout << "ERROR GetSmearMatrix: cannot find " << histoname
              << " in " << filename << "\n";
    file->Close();
    delete file;
    return nullptr;
  }

  TString name = histo->GetName();
  gROOT->cd();
  TH2F* copy = dynamic_cast<TH2F*>(histo->Clone(name));
  copy->SetDirectory(nullptr);
  file->Close();
  delete file;
  return copy;
}

DLM_Histo<double>* MakeDlmFeedDownFractionFromLambdaTH1(TH1* hLambda,
                                                        const double xScale)
{
  if(!hLambda) return nullptr;

  const int nBins = hLambda->GetNbinsX();
  double axisMin = hLambda->GetXaxis()->GetBinLowEdge(1) * xScale;
  double axisMax = hLambda->GetXaxis()->GetBinUpEdge(nBins) * xScale;

  DLM_Histo<double>* hFeedDown = new DLM_Histo<double>();
  hFeedDown->SetUp(1);
  hFeedDown->SetUp(0, nBins, axisMin, axisMax);
  hFeedDown->Initialize();

  const double weakScale = (1.0 - fracWeakDecay) * (1.0 - fracWeakDecay);
  unsigned whichBin[1];
  for(int iBin=0; iBin<nBins; ++iBin){
    const double lambdaMain = hLambda->GetBinContent(iBin+1) * weakScale;
    const double lambdaErr = hLambda->GetBinError(iBin+1) * weakScale;
    if(!std::isfinite(lambdaMain) || lambdaMain < 0.0 || lambdaMain > 1.0){
      std::cout << "ERROR: invalid lambda in bin " << iBin+1
                << ": " << lambdaMain << "\n";
      delete hFeedDown;
      return nullptr;
    }
    whichBin[0] = iBin;
    hFeedDown->SetBinContent(whichBin, 1.0 - lambdaMain);
    hFeedDown->SetBinError(whichBin, lambdaErr);
  }

  return hFeedDown;
}

DLM_Histo<double>* LoadDlmFeedDownFraction(TString fileName,
                                           TString histName,
                                           const double xScale=1000.0)
{
  TFile* file = TFile::Open(fileName, "READ");
  if(!file || file->IsZombie()){
    std::cout << "ERROR: cannot open lambda ROOT file: " << fileName << "\n";
    if(file){ file->Close(); delete file; }
    return nullptr;
  }

  TH1* hLambda = dynamic_cast<TH1*>(file->Get(histName));
  if(!hLambda){
    std::cout << "ERROR: cannot find lambda histogram: " << histName << "\n";
    file->Close();
    delete file;
    return nullptr;
  }

  DLM_Histo<double>* hFeedDown =
      MakeDlmFeedDownFractionFromLambdaTH1(hLambda, xScale);
  file->Close();
  delete file;
  return hFeedDown;
}

void CATS_pipi_Basic(CATS& kitty,
                     const unsigned numMomBins,
                     const double kMin,
                     const double kMax,
                     const bool identical=true,
                     const bool withCoulomb=true)
{
  kitty.SetMomBins(numMomBins, kMin, kMax);

  kitty.SetNumChannels(1);
  kitty.SetNumPW(0, 0);
  kitty.SetSpin(0, 0);
  kitty.SetChannelWeight(0, 1.0);

  if(withCoulomb) kitty.SetQ1Q2(1);
  else            kitty.SetQ1Q2(0);

  if(identical) kitty.SetPdgId(211, 211);
  else          kitty.SetPdgId(-211, 211);

  if(withCoulomb)    kitty.SetRedMass(0.5*Mass_Pic);
  else if(identical) kitty.SetRedMass(0.5*Mass_Pi0);
  else               kitty.SetRedMass((Mass_Pi0*Mass_Pic)/(Mass_Pi0+Mass_Pic));
}

void CATS_GaussSource(CATS& kitty, const double sourceSize)
{
  CATSparameters pars(CATSparameters::tSource, 1, true);
  pars.SetParameter(0, sourceSize);
  kitty.SetAnaSource(GaussSource, pars);
  kitty.SetUseAnalyticSource(true);
  kitty.SetAutoNormSource(false);
}

DLM_CleverMcLevyResoTM* BuildCleverMcLevyResoTMSource(CATS& kitty,
                                                      const double sourceSize,
                                                      const double alpha,
                                                      const bool fixAlpha,
                                                      const double fracReso,
                                                      const double tauReso,
                                                      const double massReso,
                                                      const double cutoff)
{
  (void)massReso;
  DLM_CleverMcLevyResoTM* source = new DLM_CleverMcLevyResoTM();

  if(fixAlpha) source->InitStability(1, alpha-1e-6, alpha+1e-6);
  else         source->InitStability(21, 1.0, 2.0);

  source->InitScale(160, 4.0, 30.0);
  source->InitRad(400, 1e-4, 128.0);
  source->InitType(2);
  source->SetUpReso(0, fracReso);
  source->SetUpReso(1, fracReso);

  Float_t k_D = 0;
  Float_t fP1 = 0;
  Float_t fP2 = 0;
  Float_t fM1 = 0;
  Float_t fM2 = 0;
  Float_t Tau1 = 0;
  Float_t Tau2 = 0;
  Float_t AngleRcP1 = 0;
  Float_t AngleRcP2 = 0;
  Float_t AngleP1P2 = 0;

  DLM_Random ranGen(11);

  TFile* filePR = TFile::Open(PipiIO::LegacyPRFileName, "READ");
  if(!filePR || filePR->IsZombie()){
    std::cout << "ERROR: cannot open PR file " << PipiIO::LegacyPRFileName << "\n";
    if(filePR){ filePR->Close(); delete filePR; }
  }
  else{
    TNtuple* tuplePR = dynamic_cast<TNtuple*>(filePR->Get("InfoTuple_ClosePairs"));
    if(!tuplePR){
      std::cout << "ERROR: cannot find InfoTuple_ClosePairs in "
                << PipiIO::LegacyPRFileName << "\n";
    }
    else{
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
      for(Long64_t iEntry=0; iEntry<tuplePR->GetEntries(); ++iEntry){
        tuplePR->GetEntry(iEntry);
        if(k_D > cutoff) continue;
        if(fP2 <= 0.f || fM2 <= 0.f || tauReso <= 0.0) continue;

        const double shift = ranGen.Exponential(fM2/(fP2*tauReso));
        source->AddBGT_PR(shift, -std::cos(AngleRcP2));
        source->AddBGT_RP(shift,  std::cos(AngleRcP2));
        ++added;
      }
      std::cout << "Loaded " << added << " PR/RP resonance kernels from "
                << PipiIO::LegacyPRFileName << "\n";
    }
    filePR->Close();
    delete filePR;
  }

  TFile* fileRR = TFile::Open(PipiIO::LegacyRRFileName, "READ");
  if(!fileRR || fileRR->IsZombie()){
    std::cout << "ERROR: cannot open RR file " << PipiIO::LegacyRRFileName << "\n";
    if(fileRR){ fileRR->Close(); delete fileRR; }
  }
  else{
    TNtuple* tupleRR = dynamic_cast<TNtuple*>(fileRR->Get("InfoTuple_ClosePairs"));
    if(!tupleRR){
      std::cout << "ERROR: cannot find InfoTuple_ClosePairs in "
                << PipiIO::LegacyRRFileName << "\n";
    }
    else{
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
      for(Long64_t iEntry=0; iEntry<tupleRR->GetEntries(); ++iEntry){
        tupleRR->GetEntry(iEntry);
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
      std::cout << "Loaded " << added << " RR resonance kernels from "
                << PipiIO::LegacyRRFileName << "\n";
    }
    fileRR->Close();
    delete fileRR;
  }

  if(fixAlpha && (std::fabs(alpha-1.0)<1e-12 || std::fabs(alpha-2.0)<1e-12)){
    source->InitNumMcIter(1000000);
  }
  else{
    source->InitNumMcIter(100000);
  }

  kitty.SetAnaSource(CatsSourceForwarder, source, 2);
  kitty.SetAnaSource(0, sourceSize);
  kitty.SetAnaSource(1, alpha);
  kitty.SetUseAnalyticSource(true);
  kitty.SetAutoNormSource(false);

  return source;
}

double femtoFuncCleverMcLevyResoTM(double* x, double* par)
{
  const double mom = x[0];
  const double sourceR = par[0];
  const double alpha = par[1];

  static double lastR = 1e99;
  static double lastAlpha = 1e99;
  static DLM_CkDecomposition* lastFitter = nullptr;
  static bool first = true;

  if(first || gFemto_pipi_fitter != lastFitter ||
     std::fabs(sourceR-lastR) > 1e-6 ||
     std::fabs(alpha-lastAlpha) > 1e-6){
    lastR = sourceR;
    lastAlpha = alpha;
    lastFitter = gFemto_pipi_fitter;
    first = false;

    gFemto_pipi_fitter->GetCk()->SetSourcePar(0, sourceR);
    gFemto_pipi_fitter->GetCk()->SetSourcePar(1, alpha);
    gFemto_pipi_fitter->Update(true);
  }

  return gFemto_pipi_fitter->EvalCk(mom);
}

double femtoFuncGaussOnly(double* x, double* par)
{
  const double mom = x[0];
  const double sourceR = par[0];

  static double lastR = 1e99;
  static DLM_CkDecomposition* lastFitter = nullptr;
  static bool first = true;

  if(first || gFemto_pipi_fitter != lastFitter ||
     std::fabs(sourceR-lastR) > 1e-6){
    lastR = sourceR;
    lastFitter = gFemto_pipi_fitter;
    first = false;

    gFemto_pipi_fitter->GetCk()->SetSourcePar(0, sourceR);
    gFemto_pipi_fitter->Update(true);
  }

  return gFemto_pipi_fitter->EvalCk(mom);
}

double FitFuncCleverMcLevyResoTMPol1(double* x, double* par)
{
  if(!gKittyPiPi || !gCkPiPi || !gFemto_pipi_fitter) return 1.0;

  const double mom = x[0];
  const double ck = femtoFuncCleverMcLevyResoTM(x, &par[0]);
  const double baseline = par[2] * (1.0 + par[3]*mom);
  return ck * baseline;
}

double FitFuncGaussOnlyPol1(double* x, double* par)
{
  if(!gKittyPiPi || !gCkPiPi || !gFemto_pipi_fitter) return 1.0;

  const double mom = x[0];
  const double ck = femtoFuncGaussOnly(x, &par[0]);
  const double baseline = par[1] * (1.0 + par[2]*mom);
  return ck * baseline;
}

double baseline3(double* x, double* p)
{
  const double xx = x[0];
  return p[0] * (1.0 + p[1]*xx + p[2]*xx*xx + p[3]*xx*xx*xx);
}

void DrawErrorBand(TH1* hData,
                   TF1* fit,
                   TGraphErrors*& grfit,
                   TGraphAsymmErrors*& grBand,
                   TGraph& gUp,
                   TGraph& gDn,
                   TGraph& gCtrl,
                   double CL = 0.95)
{
  (void)grfit;
  (void)grBand;
  (void)gUp;
  (void)gDn;
  (void)gCtrl;
  (void)CL;

  const double xmin = fit->GetXmin();
  const double xmax = fit->GetXmax();
  TFitResultPtr r = hData->Fit(fit, "RS+", "", xmin, xmax);
  if(!r.Get() || !r->IsValid()){
    std::cout << "WARNING: fit result is not valid.\n";
  }

  fit->SetLineColor(kRed);
  fit->SetLineStyle(2);
  fit->SetLineWidth(2);
}

void BuildFitSourceRadiusBand(const TH1* hTemplate,
                              TF1* fit,
                              const int sourcePar,
                              const double sourceErr,
                              const TString& prefix,
                              TH1*& hCentral,
                              TH1*& hUp,
                              TH1*& hDown,
                              TGraphAsymmErrors*& grBand)
{
  hCentral = nullptr;
  hUp = nullptr;
  hDown = nullptr;
  grBand = nullptr;
  if(!hTemplate || !fit || sourcePar < 0) return;

  const double bestR = fit->GetParameter(sourcePar);
  const double errR = std::fabs(sourceErr);
  const double loR = std::max(1e-6, bestR - errR);
  const double hiR = bestR + errR;

  hCentral = dynamic_cast<TH1*>(hTemplate->Clone(prefix + "_central"));
  hUp = dynamic_cast<TH1*>(hTemplate->Clone(prefix + "_sourceR_plus1sigma"));
  hDown = dynamic_cast<TH1*>(hTemplate->Clone(prefix + "_sourceR_minus1sigma"));
  if(!hCentral || !hUp || !hDown){
    delete hCentral;
    delete hUp;
    delete hDown;
    hCentral = nullptr;
    hUp = nullptr;
    hDown = nullptr;
    return;
  }

  hCentral->SetDirectory(nullptr);
  hUp->SetDirectory(nullptr);
  hDown->SetDirectory(nullptr);
  hCentral->Reset("ICES");
  hUp->Reset("ICES");
  hDown->Reset("ICES");

  grBand = new TGraphAsymmErrors(hTemplate->GetNbinsX());
  grBand->SetName(prefix + "_sourceR_1sigma_band");
  grBand->SetTitle(prefix + "_sourceR_1sigma_band");

  const int nBins = hTemplate->GetNbinsX();
  std::vector<double> xVals(nBins, 0.0);
  std::vector<double> exVals(nBins, 0.0);
  std::vector<double> yCentralVals(nBins, 0.0);
  std::vector<double> yHiRVals(nBins, 0.0);
  std::vector<double> yLoRVals(nBins, 0.0);

  fit->SetParameter(sourcePar, bestR);
  for(int iBin=1; iBin<=nBins; ++iBin){
    const int idx = iBin-1;
    xVals[idx] = hTemplate->GetXaxis()->GetBinCenter(iBin);
    exVals[idx] = 0.5*hTemplate->GetXaxis()->GetBinWidth(iBin);
    yCentralVals[idx] = fit->Eval(xVals[idx]);
  }

  fit->SetParameter(sourcePar, hiR);
  for(int iBin=1; iBin<=nBins; ++iBin){
    const int idx = iBin-1;
    yHiRVals[idx] = fit->Eval(xVals[idx]);
  }

  fit->SetParameter(sourcePar, loR);
  for(int iBin=1; iBin<=nBins; ++iBin){
    const int idx = iBin-1;
    yLoRVals[idx] = fit->Eval(xVals[idx]);
  }

  for(int iBin=1; iBin<=nBins; ++iBin){
    const int idx = iBin-1;
    const double x = hTemplate->GetXaxis()->GetBinCenter(iBin);
    const double ex = 0.5*hTemplate->GetXaxis()->GetBinWidth(iBin);
    const double yCentral = yCentralVals[idx];
    const double yHiR = yHiRVals[idx];
    const double yLoR = yLoRVals[idx];

    const double yMin = std::min(yHiR, yLoR);
    const double yMax = std::max(yHiR, yLoR);

    hCentral->SetBinContent(iBin, yCentral);
    hUp->SetBinContent(iBin, yMax);
    hDown->SetBinContent(iBin, yMin);
    hCentral->SetBinError(iBin, 0.0);
    hUp->SetBinError(iBin, 0.0);
    hDown->SetBinError(iBin, 0.0);

    grBand->SetPoint(iBin-1, x, yCentral);
    grBand->SetPointError(iBin-1, ex, ex,
                          std::max(0.0, yCentral-yMin),
                          std::max(0.0, yMax-yCentral));
  }

  fit->SetParameter(sourcePar, bestR);
}

} // namespace

void Init_pipi_Theory(TString SourceType,
                      const bool Identical,
                      const bool WithCoulomb,
                      const bool onlyTheory,
                      const double SourceSize,
                      const double SourceAlpha,
                      const unsigned NumMomBins,
                      const double kMin,
                      const double kMax,
                      TString lambdaFileName,
                      TString lambdaHistName)
{
  (void)onlyTheory;

  CleanupPipiTheory();
  gKittyPiPi = new CATS();

  const bool usePlainGauss = SourceType == "Gauss";
  const bool fixAlpha =
      SourceType != "CoreResoFreeAlpha" &&
      SourceType != "DLM_CleverMcLevyResoTM_FreeAlpha";
  if(SourceType!="Gauss" &&
     SourceType!="CoreReso" &&
     SourceType!="CoreResoFixedAlpha" &&
     SourceType!="CoreResoFreeAlpha" &&
     SourceType!="DLM_CleverMcLevyResoTM" &&
     SourceType!="DLM_CleverMcLevyResoTM_FreeAlpha"){
    std::cout << "WARNING Init_pipi_Theory: source '" << SourceType
              << "' is disabled in the cleaned ExtendedCk. "
              << "Using DLM_CleverMcLevyResoTM instead.\n";
  }

  if(usePlainGauss){
    CATS_GaussSource(*gKittyPiPi, SourceSize);
  }
  else{
    gCleverMcLevyResoTM_pipi = BuildCleverMcLevyResoTMSource(
        *gKittyPiPi,
        SourceSize,
        SourceAlpha,
        fixAlpha,
        fracShortLivedResonance,
        1.5,
        1124.0,
        kMax);
  }

  CATS_pipi_Basic(*gKittyPiPi, NumMomBins, kMin, kMax, Identical, WithCoulomb);

  gKittyPiPi->SetNotifications(CATS::nWarning);
  gKittyPiPi->SetNotifications(CATS::nError);
  gKittyPiPi->KillTheCat();

  gCkPiPi = new DLM_Ck(gKittyPiPi->GetNumSourcePars(), 0, *gKittyPiPi);
  gCkPiPi->Update(true);

  TH2F* hMomentumResolution_pipi =
      GetSmearMatrix(PipiIO::SmearingFileName, PipiIO::SmearingHistoName);
  if(!hMomentumResolution_pipi){
    std::cout << "ERROR: failed to get momentum smearing matrix\n";
    return;
  }

  DLM_Histo<double>* hLambdaFeedDown =
      LoadDlmFeedDownFraction(lambdaFileName, lambdaHistName, 1000.0);
  if(!hLambdaFeedDown){
    std::cout << "ERROR: failed to build k*-dependent lambda for AddContribution\n";
    delete hMomentumResolution_pipi;
    return;
  }

  DLM_CkDecomposition* ckDec =
      new DLM_CkDecomposition("pipi", 1, *gCkPiPi, hMomentumResolution_pipi);
  ckDec->AddContribution(0, *hLambdaFeedDown, DLM_CkDecomposition::cFeedDown);

  delete hLambdaFeedDown;

  ckDec->Update(true, true);
  gFemto_pipi_fitter = ckDec;
}

void fitDataCoreReso(bool Identical,
                     bool WithCoulomb,
                     double alphaInit,
                     bool fixAlpha)
{
  TFile* input = TFile::Open(PipiIO::DataFileName, "READ");
  if(!input || input->IsZombie()){
    std::cout << "ERROR: cannot open input ROOT file: "
              << PipiIO::DataFileName << "\n";
    if(input){ input->Close(); delete input; }
    return;
  }

  TH1* hData = dynamic_cast<TH1*>(input->Get(PipiIO::DataHistoName));
  if(!hData){
    std::cout << "ERROR: cannot find histogram: "
              << PipiIO::DataHistoName << "\n";
    input->Close();
    delete input;
    return;
  }

  TH1* hDataToFit = dynamic_cast<TH1*>(hData->Clone("copyhisto_CleverMcLevyResoTM"));
  if(!hDataToFit){
    std::cout << "ERROR: cannot clone input histogram\n";
    input->Close();
    delete input;
    return;
  }
  hDataToFit->SetDirectory(nullptr);

  const double sourceRInit = 12.0;
  const TString sourceType =
      fixAlpha ? "DLM_CleverMcLevyResoTM" : "DLM_CleverMcLevyResoTM_FreeAlpha";
  Init_pipi_Theory(sourceType,
                   Identical,
                   WithCoulomb,
                   false,
                   sourceRInit,
                   alphaInit,
                   PipiIO::CatsNumMomBins,
                   PipiIO::CatsKMin,
                   PipiIO::CatsKMax,
                   PipiIO::LambdaFileName,
                   PipiIO::LambdaHistoName);
  if(!gFemto_pipi_fitter){
    std::cout << "ERROR: Init_pipi_Theory failed, skip fit.\n";
    input->Close();
    delete input;
    return;
  }

  TF1* baselinePrefit =
      new TF1("baseline_prefit_CleverMcLevyResoTM",
              "[0]*(1.0+[1]*x)",
              PipiIO::BaselineFitMin,
              PipiIO::BaselineFitMax);
  baselinePrefit->SetParameter(0, 1.0);
  baselinePrefit->SetParameter(1, 0.0);
  hDataToFit->Fit(baselinePrefit, "RQ0", "",
                  PipiIO::BaselineFitMin,
                  PipiIO::BaselineFitMax);

  TF1* fit = new TF1("fit_pipi_CleverMcLevyResoTM",
                     FitFuncCleverMcLevyResoTMPol1,
                     PipiIO::FitMin,
                     PipiIO::FitMax,
                     4);
  fit->SetNpx(PipiIO::FitNpx);
  fit->SetParName(0, "R_{Levy}");
  fit->SetParName(1, "#alpha_{Levy}");
  fit->SetParName(2, "baseline p_{0}");
  fit->SetParName(3, "baseline p_{1}");

  fit->SetParameter(0, sourceRInit);
  fit->SetParLimits(0, 3.0, 25.0);

  fit->SetParameter(1, alphaInit);
  if(fixAlpha){
    fit->FixParameter(1, alphaInit);
  }
  else{
    fit->SetParLimits(1, 1.0, 2.0);
  }

  fit->SetParameter(2, baselinePrefit->GetParameter(0));
  fit->FixParameter(2, baselinePrefit->GetParameter(0));
  fit->SetParameter(3, baselinePrefit->GetParameter(1));
  fit->FixParameter(3, baselinePrefit->GetParameter(1));

  TGraphAsymmErrors* grBand = new TGraphAsymmErrors();
  TGraphErrors* grfit = new TGraphErrors();
  TGraph gUp;
  TGraph gDn;
  TGraph gCtrl;

  std::cout << "start DLM_CleverMcLevyResoTM fit: alpha="
            << alphaInit << (fixAlpha ? " fixed" : " free") << "\n";
  DrawErrorBand(hDataToFit, fit, grfit, grBand, gUp, gDn, gCtrl, 0.95);

  std::cout << "DLM_CleverMcLevyResoTM fit done:\n"
            << "  R = " << fit->GetParameter(0)
            << " +/- " << fit->GetParError(0) << "\n"
            << "  alpha = " << fit->GetParameter(1)
            << " +/- " << fit->GetParError(1)
            << (fixAlpha ? " (fixed)" : "") << "\n"
            << "  chi2/ndf = " << fit->GetChisquare() << "/"
            << fit->GetNDF();
  if(fit->GetNDF()>0){
    std::cout << " = " << fit->GetChisquare()/double(fit->GetNDF());
  }
  std::cout << "\n";

  TFile* output = new TFile(PipiIO::OutputFitCoreReso1D, "RECREATE");
  output->cd();
  hDataToFit->Write();
  baselinePrefit->Write();

  fit->SetName("CATS_fit_CleverMcLevyResoTM");
  fit->Write();
  // Backward-compatible name for old result readers.
  fit->Write("CATS_fit_coreReso", TObject::kOverwrite);

  TH1* hFitCentral = nullptr;
  TH1* hFitUp = nullptr;
  TH1* hFitDown = nullptr;
  TGraphAsymmErrors* grSourceBand = nullptr;
  BuildFitSourceRadiusBand(hDataToFit,
                           fit,
                           0,
                           fit->GetParError(0),
                           "CATS_fit_CleverMcLevyResoTM",
                           hFitCentral,
                           hFitUp,
                           hFitDown,
                           grSourceBand);
  if(hFitCentral) hFitCentral->Write();
  if(hFitUp) hFitUp->Write();
  if(hFitDown) hFitDown->Write();
  if(grSourceBand) grSourceBand->Write();

  TF1* baseline = new TF1("fit_pipi_baseline_CleverMcLevyResoTM",
                          baseline3, 0, 400, 4);
  baseline->SetParameter(0, fit->GetParameter(2));
  baseline->SetParameter(1, fit->GetParameter(3));
  baseline->SetParameter(2, 0.0);
  baseline->SetParameter(3, 0.0);
  baseline->SetLineColor(kGreen-2);
  baseline->Write();

  output->Close();
  input->Close();

  std::cout << "Saved DLM_CleverMcLevyResoTM fit: "
            << PipiIO::OutputFitCoreReso1D << "\n";

  delete grBand;
  delete grfit;
  delete grSourceBand;
  delete hFitCentral;
  delete hFitUp;
  delete hFitDown;
  delete baseline;
  delete baselinePrefit;
  delete fit;
  delete output;
  delete hDataToFit;
  delete input;
}

void fitDataCoreResoFixedAlpha(bool Identical,
                               bool WithCoulomb,
                               double alpha)
{
  const TString savedOutput = PipiIO::OutputFitCoreReso1D;
  TString alphaTag = Form("%.1f", alpha);
  alphaTag.ReplaceAll(".", "p");
  PipiIO::OutputFitCoreReso1D =
      Form("./out/Aug/fit_results_CleverMcLevyResoTM_ForMax_fixedAlpha%s_cent10to30_mT4.root",
           alphaTag.Data());
  fitDataCoreReso(Identical, WithCoulomb, alpha, true);
  PipiIO::OutputFitCoreReso1D = savedOutput;
}

void fitDataCoreResoFreeAlpha(bool Identical,
                              bool WithCoulomb,
                              double alphaInit)
{
  fitDataCoreReso(Identical, WithCoulomb, alphaInit, false);
}

void fitData(bool Identical, bool WithCoulomb)
{
  std::cout << "fitData alias: running DLM_CleverMcLevyResoTM with free alpha.\n";
  fitDataCoreResoFreeAlpha(Identical, WithCoulomb, sourceAlpha);
}

void fitDataGaussOnly(bool Identical, bool WithCoulomb)
{
  TFile* input = TFile::Open(PipiIO::DataFileName, "READ");
  if(!input || input->IsZombie()){
    std::cout << "ERROR: cannot open input ROOT file: "
              << PipiIO::DataFileName << "\n";
    if(input){ input->Close(); delete input; }
    return;
  }

  TH1* hData = dynamic_cast<TH1*>(input->Get(PipiIO::DataHistoName));
  if(!hData){
    std::cout << "ERROR: cannot find histogram: "
              << PipiIO::DataHistoName << "\n";
    input->Close();
    delete input;
    return;
  }

  TH1* hDataToFit = dynamic_cast<TH1*>(hData->Clone("copyhisto_gaussOnly"));
  if(!hDataToFit){
    std::cout << "ERROR: cannot clone input histogram\n";
    input->Close();
    delete input;
    return;
  }
  hDataToFit->SetDirectory(nullptr);

  const double sourceRInit = 12.0;
  Init_pipi_Theory("Gauss",
                   Identical,
                   WithCoulomb,
                   false,
                   sourceRInit,
                   2.0,
                   PipiIO::CatsNumMomBins,
                   PipiIO::CatsKMin,
                   PipiIO::CatsKMax,
                   PipiIO::LambdaFileName,
                   PipiIO::LambdaHistoName);
  if(!gFemto_pipi_fitter){
    std::cout << "ERROR: Init_pipi_Theory failed, skip Gauss-only fit.\n";
    input->Close();
    delete input;
    return;
  }

  TF1* baselinePrefit =
      new TF1("baseline_prefit_gaussOnly",
              "[0]*(1.0+[1]*x)",
              PipiIO::BaselineFitMin,
              PipiIO::BaselineFitMax);
  baselinePrefit->SetParameter(0, 1.0);
  baselinePrefit->SetParameter(1, 0.0);
  hDataToFit->Fit(baselinePrefit, "RQ0", "",
                  PipiIO::BaselineFitMin,
                  PipiIO::BaselineFitMax);

  TF1* fit = new TF1("fit_pipi_gaussOnly",
                     FitFuncGaussOnlyPol1,
                     PipiIO::FitMin,
                     PipiIO::FitMax,
                     3);
  fit->SetNpx(PipiIO::FitNpx);
  fit->SetParName(0, "R_{Gauss}");
  fit->SetParName(1, "baseline p_{0}");
  fit->SetParName(2, "baseline p_{1}");
  fit->SetParameter(0, sourceRInit);
  fit->SetParLimits(0, 3.0, 25.0);
  fit->SetParameter(1, baselinePrefit->GetParameter(0));
  fit->FixParameter(1, baselinePrefit->GetParameter(0));
  fit->SetParameter(2, baselinePrefit->GetParameter(1));
  fit->FixParameter(2, baselinePrefit->GetParameter(1));

  TGraphAsymmErrors* grBand = new TGraphAsymmErrors();
  TGraphErrors* grfit = new TGraphErrors();
  TGraph gUp;
  TGraph gDn;
  TGraph gCtrl;

  std::cout << "start Gauss-only fit\n";
  DrawErrorBand(hDataToFit, fit, grfit, grBand, gUp, gDn, gCtrl, 0.95);

  std::cout << "Gauss-only fit done:\n"
            << "  R = " << fit->GetParameter(0)
            << " +/- " << fit->GetParError(0) << "\n"
            << "  chi2/ndf = " << fit->GetChisquare() << "/"
            << fit->GetNDF();
  if(fit->GetNDF()>0){
    std::cout << " = " << fit->GetChisquare()/double(fit->GetNDF());
  }
  std::cout << "\n";

  TFile* output = new TFile(PipiIO::OutputFitGaussOnly1D, "RECREATE");
  output->cd();
  hDataToFit->Write();
  baselinePrefit->Write();
  fit->SetName("CATS_fit_gaussOnly");
  fit->Write();

  TH1* hFitCentral = nullptr;
  TH1* hFitUp = nullptr;
  TH1* hFitDown = nullptr;
  TGraphAsymmErrors* grSourceBand = nullptr;
  BuildFitSourceRadiusBand(hDataToFit,
                           fit,
                           0,
                           fit->GetParError(0),
                           "CATS_fit_gaussOnly",
                           hFitCentral,
                           hFitUp,
                           hFitDown,
                           grSourceBand);
  if(hFitCentral) hFitCentral->Write();
  if(hFitUp) hFitUp->Write();
  if(hFitDown) hFitDown->Write();
  if(grSourceBand) grSourceBand->Write();

  TF1* baseline = new TF1("fit_pipi_baseline_gaussOnly",
                          baseline3, 0, 400, 4);
  baseline->SetParameter(0, fit->GetParameter(1));
  baseline->SetParameter(1, fit->GetParameter(2));
  baseline->SetParameter(2, 0.0);
  baseline->SetParameter(3, 0.0);
  baseline->SetLineColor(kGreen-2);
  baseline->Write();

  output->Close();
  input->Close();

  std::cout << "Saved Gauss-only fit: "
            << PipiIO::OutputFitGaussOnly1D << "\n";

  delete grBand;
  delete grfit;
  delete grSourceBand;
  delete hFitCentral;
  delete hFitUp;
  delete hFitDown;
  delete baseline;
  delete baselinePrefit;
  delete fit;
  delete output;
  delete hDataToFit;
  delete input;
}

void fitData3DTo1D(bool Identical, bool WithCoulomb, bool usePRF3DShift)
{
  (void)usePRF3DShift;
  std::cout << "fitData3DTo1D alias: 3D-to-1D source is disabled; "
            << "running DLM_CleverMcLevyResoTM with fixed alpha.\n";
  fitDataCoreResoFixedAlpha(Identical, WithCoulomb, sourceAlpha);
}

void Ck_pL_Ledni_Usmani()
{
  std::cout << "Ck_pL_Ledni_Usmani is not part of the cleaned pion fit code.\n";
}

void RootFile_DlmCk(const TString& RootFileName,
                    const TString& GraphName,
                    DLM_Ck* CkToPlot)
{
  if(!CkToPlot){
    std::cout << "ERROR RootFile_DlmCk: null DLM_Ck\n";
    return;
  }

  TFile file(RootFileName, "RECREATE");
  const unsigned nBins = CkToPlot->GetNbins();
  TGraph graph(nBins);
  graph.SetName(GraphName);
  for(unsigned iBin=0; iBin<nBins; ++iBin){
    graph.SetPoint(iBin,
                   CkToPlot->GetBinCenter(0, iBin)/1000.0,
                   CkToPlot->GetBinContent(iBin));
  }
  graph.Write();
}

void RootFile_DlmCk(const TString& RootFileName,
                    const TString& GraphName,
                    DLM_CkDecomposition* CkToPlot,
                    bool PlotIndividualContributions)
{
  (void)PlotIndividualContributions;
  if(!CkToPlot){
    std::cout << "ERROR RootFile_DlmCk: null DLM_CkDecomposition\n";
    return;
  }

  TFile file(RootFileName, "RECREATE");
  const unsigned nBins = CkToPlot->GetCk()->GetNbins();
  TGraph graph(nBins);
  graph.SetName(GraphName);
  for(unsigned iBin=0; iBin<nBins; ++iBin){
    const double k = CkToPlot->GetCk()->GetBinCenter(0, iBin);
    graph.SetPoint(iBin, k/1000.0, CkToPlot->EvalCk(k));
  }
  graph.Write();
}

void Ck_pp_Decomposition(const TString& SourceType)
{
  std::cout << "Ck_pp_Decomposition is disabled in the cleaned pion fit code. "
            << "Requested source: " << SourceType << "\n";
}
