#include "Basics.h"
#include "ExtendedCk.h"

#include "/Users/wuwenya/Documents/code/femto_flow/CATS/CATS_install/install/include/CATS.h"
#include "/Users/wuwenya/Documents/code/femto_flow/CATS/CATS_install/install/include/CATStools.h"
#include "/Users/wuwenya/Documents/code/femto_flow/CATS/CATS_install/install/include/CATSconstants.h"

#include "/Users/wuwenya/Documents/code/femto_flow/CATS/CATS_install/install/include/DLM_Source.h"
#include "/Users/wuwenya/Documents/code/femto_flow/CATS/CATS_install/install/include/DLM_Ck.h"
#include "/Users/wuwenya/Documents/code/femto_flow/CATS/CATS_install/install/include/DLM_CkDecomposition.h"
#include "/Users/wuwenya/Documents/code/femto_flow/CATS/CATS_install/install/include/DLM_Histo.h"

#include "TFile.h"
#include "TH1.h"
#include "TH1D.h"
#include "TF1.h"
#include "TH2F.h"
#include "TTree.h"
#include "TString.h"
#include "TROOT.h"
#include "TNtuple.h"
#include "TGraph.h"
#include "TRandom3.h"
#include "TMath.h"
#include "TGraphAsymmErrors.h"
#include "TGraphErrors.h"
#include "TFitResultPtr.h"
#include "TFitResult.h"

#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>


enum class FlowMode {
  None,
  V2YieldOnly,
  V2YieldAndGeometry
};


namespace PipiIO {
  const TString DataFileName = "../Files/cf/cf_data_cent0to10_tofOnly_PbPb.root";
  const TString DataHistoName = "cf_same_proj_mt0_cent0_scaled";

  const TString LambdaFileName = "../Files/lam/kStarDepLambda_cent0to10_pion_PbPb.root";
  const TString LambdaHistoName = "hLambdaKstar_mT0.30to0.50_cent0to10";

  const TString SmearingFileName = "../Files/resolution_rebin4.root";
  const TString SmearingHistoName = "resolution";

  const TString LegacyPRFileName = "../Files/ForMax_pi_piReso_noOmega.root";
  const TString LegacyRRFileName = "../Files/ForMax_piReso_piReso_noOmega.root";
  const TString PRF3DPRFileName = "../Files/epos3D_pi_piReso_noOmega.root";
  const TString PRF3DRRFileName = "../Files/epos3D_piReso_piReso_noOmega.root";

  const TString OutputFit1D = "./out/fit_results_flxAlpha_1D_cent0to10.root";
  const TString OutputFit3DTo1D = "./out/fit_results_flxAlpha_3DTo1D_cent0to10.root";
  const TString OutputScanAlphaBinning = "scan_alpha_binning_fixedWeak.root";
  const TString OutputScanAlpha8MeV = "scan_alpha8MeV_fixedWeak.root";
  const TString OutputScanWeakDecay8MeV = "scan_weakDecay8MeV.root";

  const unsigned CatsNumMomBins = 500;
  const double CatsKMin = 0.0;
  const double CatsKMax = 400.0;
  const int FitNpx = 500;
  const double FitMin = 2.0;
  const double FitMax = 200.0;
  const double BaselineFitMin = 200.0;
  const double BaselineFitMax = 400.0;
}

// =========================================================
// Globals (must live during the whole fit)
// =========================================================
static CATS*   gKittyPiPi = nullptr;
static DLM_Ck* gCkPiPi    = nullptr;
static DLM_CleverMcLevyResoTM* gCRSpipi = nullptr;
static DLM_CleverMcLevyResoTM_3DTo1D* gCRSpipi3DTo1D = nullptr;
static DLM_CleverMcDoubleGaussResoTM* gCRSpipiDoubleGaus = nullptr;
static DLM_CkDecomposition* gFemto_pipi_fitter = nullptr;
double fracWeakDecay = 0.12;
const double fracShortLivedResonance = 0.72;
const double doubleGaussTailWidth = 20.0;
const double doubleGaussCoreFraction = 0.80;
double sourceAlpha = 2.0;
static bool gUseDoubleGaussianSource = false;

void DrawErrorBand(TH1* hData, TF1* fit,
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
  if (!r.Get() || !r->IsValid()) {
    printf("Invalid fit result.\n");
    return;
  }

  fit->SetLineColor(kRed);
  fit->SetLineStyle(2);
  fit->SetLineWidth(2);
}

TH2F* GetSmearMatrix(TString filename, TString histoname){
    TFile* FileROOT = new TFile(filename, "read");
    TH2F* histo = (TH2F*)FileROOT->Get(histoname);
    if(!histo){printf("\033[1;31mERROR:\033[0m The histo '%s' in file '%s' does not exist\n",histoname.Data(),filename.Data());return NULL;}
    TString Name = histo->GetName();
    gROOT->cd();
    TH2F *histoCopy = (TH2F*)histo->Clone("histoCopy");
    delete FileROOT;
    histoCopy->SetName(Name);
    return histoCopy;
}

// forward declaration (important for TF1 template)
double MickeyFitter_pipi_Pol3(double* x, double* par);

// save a DLM_Ck to a file in the form of a TH1D
// NOTE: assumes uniform binning in the DLM_Ck object (true for standard CATS usage)
void RootFile_DlmCk_Histo(const TString& RootFileName,
                          const TString& HistoName,
                          DLM_Ck* CkToPlot)
{
    if(!CkToPlot){
        printf("ERROR RootFile_DlmCk_Histo: CkToPlot is null\n");
        return;
    }

    TFile* RootFile = TFile::Open(RootFileName, "UPDATE");
    if(!RootFile || RootFile->IsZombie()){
        if(RootFile){ RootFile->Close(); delete RootFile; }
        RootFile = TFile::Open(RootFileName, "RECREATE");
        if(!RootFile || RootFile->IsZombie()){
            printf("ERROR RootFile_DlmCk_Histo: cannot open output file %s\n", RootFileName.Data());
            if(RootFile){ RootFile->Close(); delete RootFile; }
            return;
        }
    }

    const unsigned NumBins = CkToPlot->GetNbins();
    if(NumBins<1){
        printf("ERROR RootFile_DlmCk_Histo: NumBins < 1\n");
        RootFile->Close(); delete RootFile;
        return;
    }

    // infer bin width from first two bin centers (uniform binning assumed)
    const double x0 = CkToPlot->GetBinCenter(0,0);
    double bw = 1.0;
    if(NumBins>1){
        const double x1 = CkToPlot->GetBinCenter(0,1);
        bw = x1 - x0;
    }

    const double xMin = x0 - 0.5*bw;
    const double xMax = xMin + NumBins*bw;

    TH1D h(HistoName, HistoName, (int)NumBins, xMin/1000., xMax/1000.);
    h.Sumw2(false); // no errors provided by DLM_Ck; keep simple

    for(unsigned uBin=0; uBin<NumBins; uBin++){
        h.SetBinContent((int)uBin+1, CkToPlot->GetBinContent(uBin));
    }

    RootFile->cd();
    h.Write("", TObject::kOverwrite);

    RootFile->Close();
    delete RootFile;
}

template <class SourceType>
void LoadOneDimensionalResonanceKernels(SourceType* source,
                                        const double resonanceTau,
                                        const double cutoff)
{
  Float_t k_D = 0;
  Float_t fP1 = 0;
  Float_t fP2 = 0;
  Float_t fM1 = 0;
  Float_t fM2 = 0;
  Float_t tau1 = 0;
  Float_t tau2 = 0;
  Float_t angleRcP1 = 0;
  Float_t angleRcP2 = 0;
  Float_t angleP1P2 = 0;

  DLM_Random ranGen(11);

  TFile* filePR = TFile::Open(PipiIO::LegacyPRFileName, "READ");
  if(!filePR || filePR->IsZombie()) {
    std::cout << "ERROR: cannot open primary-resonance file: "
              << PipiIO::LegacyPRFileName << "\n";
    if(filePR) {
      filePR->Close();
      delete filePR;
    }
  }
  else {
    TNtuple* tuplePR = dynamic_cast<TNtuple*>(filePR->Get("InfoTuple_ClosePairs"));
    if(!tuplePR) {
      std::cout << "ERROR: cannot find InfoTuple_ClosePairs in "
                << PipiIO::LegacyPRFileName << "\n";
    }
    else {
      tuplePR->SetBranchAddress("k_D", &k_D);
      tuplePR->SetBranchAddress("P1", &fP1);
      tuplePR->SetBranchAddress("P2", &fP2);
      tuplePR->SetBranchAddress("M1", &fM1);
      tuplePR->SetBranchAddress("M2", &fM2);
      tuplePR->SetBranchAddress("Tau1", &tau1);
      tuplePR->SetBranchAddress("Tau2", &tau2);
      tuplePR->SetBranchAddress("AngleRcP1", &angleRcP1);
      tuplePR->SetBranchAddress("AngleRcP2", &angleRcP2);
      tuplePR->SetBranchAddress("AngleP1P2", &angleP1P2);

      unsigned addedPR = 0;
      for(Long64_t iEntry = 0; iEntry < tuplePR->GetEntries(); ++iEntry) {
        tuplePR->GetEntry(iEntry);
        if(k_D > cutoff) continue;
        if(fP2 <= 0.f || resonanceTau <= 0.) continue;

        const double shift = ranGen.Exponential(fM2 / (fP2 * resonanceTau));
        source->AddBGT_PR(shift, -std::cos(angleRcP2));
        source->AddBGT_RP(shift,  std::cos(angleRcP2));
        ++addedPR;
      }
      std::cout << "Loaded " << addedPR << " PR/RP resonance kernels from "
                << PipiIO::LegacyPRFileName << "\n";
    }
    filePR->Close();
    delete filePR;
  }

  TFile* fileRR = TFile::Open(PipiIO::LegacyRRFileName, "READ");
  if(!fileRR || fileRR->IsZombie()) {
    std::cout << "ERROR: cannot open resonance-resonance file: "
              << PipiIO::LegacyRRFileName << "\n";
    if(fileRR) {
      fileRR->Close();
      delete fileRR;
    }
  }
  else {
    TNtuple* tupleRR = dynamic_cast<TNtuple*>(fileRR->Get("InfoTuple_ClosePairs"));
    if(!tupleRR) {
      std::cout << "ERROR: cannot find InfoTuple_ClosePairs in "
                << PipiIO::LegacyRRFileName << "\n";
    }
    else {
      tupleRR->SetBranchAddress("k_D", &k_D);
      tupleRR->SetBranchAddress("P1", &fP1);
      tupleRR->SetBranchAddress("P2", &fP2);
      tupleRR->SetBranchAddress("M1", &fM1);
      tupleRR->SetBranchAddress("M2", &fM2);
      tupleRR->SetBranchAddress("Tau1", &tau1);
      tupleRR->SetBranchAddress("Tau2", &tau2);
      tupleRR->SetBranchAddress("AngleRcP1", &angleRcP1);
      tupleRR->SetBranchAddress("AngleRcP2", &angleRcP2);
      tupleRR->SetBranchAddress("AngleP1P2", &angleP1P2);

      unsigned addedRR = 0;
      for(Long64_t iEntry = 0; iEntry < tupleRR->GetEntries(); ++iEntry) {
        tupleRR->GetEntry(iEntry);
        if(k_D > cutoff) continue;
        if(fM1 <= 0.f || fM2 <= 0.f) continue;
        if(fP1 <= 0.f || fP2 <= 0.f || resonanceTau <= 0.) continue;

        const double shift1 = ranGen.Exponential(fM1 / (fP1 * resonanceTau));
        const double shift2 = ranGen.Exponential(fM2 / (fP2 * resonanceTau));

        source->AddBGT_RR(
            shift1, std::cos(angleRcP1),
            shift2, std::cos(angleRcP2),
            std::cos(angleP1P2)
        );
        ++addedRR;
      }
      std::cout << "Loaded " << addedRR << " RR resonance kernels from "
                << PipiIO::LegacyRRFileName << "\n";
    }
    fileRR->Close();
    delete fileRR;
  }
}

// =========================================================
// Helper: Gaussian source for CATS
// =========================================================
void CATS_GaussSource(CATS& Kitty, const double& SourceSize){
  CATSparameters cPars(CATSparameters::tSource, 1, true);
  cPars.SetParameter(0, SourceSize);
  Kitty.SetAnaSource(GaussSource, cPars);
  Kitty.SetUseAnalyticSource(true);
  Kitty.SetAutoNormSource(false);
}

DLM_CleverMcLevyResoTM* CATS_ResoSource_pipi_UsingPionAngles(
    CATS& Kitty,
    const double& SourceSize,
    const double& Alpha,
    const bool& FixAlpha,
    const double& fracreso,
    const double& taureso,   // effective c*tau in fm (e.g. ~1.5 fm for "short-lived")
    const double& massreso,  // effective resonance mass in GeV/c^2 (e.g. ~1.1 GeV)
    const double& cutoff     // k_D cutoff in MeV/c (e.g. 200)
){
    (void)massreso;

    DLM_CleverMcLevyResoTM* CleverMcLevyResoTM_pipi = new DLM_CleverMcLevyResoTM();

    if(FixAlpha) CleverMcLevyResoTM_pipi->InitStability(1, Alpha-1e-6, Alpha+1e-6);
    else         CleverMcLevyResoTM_pipi->InitStability(21, 1, 2);

    CleverMcLevyResoTM_pipi->InitScale(160, 4.0, 30.0);
    CleverMcLevyResoTM_pipi->InitRad(400, 1e-4, 128);
    CleverMcLevyResoTM_pipi->InitType(2);

    CleverMcLevyResoTM_pipi->SetUpReso(0, fracreso);
    CleverMcLevyResoTM_pipi->SetUpReso(1, fracreso);

    Float_t k_D, fP1, fP2, fM1, fM2, Tau1, Tau2;
    Float_t AngleRcP1, AngleRcP2, AngleP1P2;

    double v2_cent = 0.035;
    TRandom3 rnd1;
    TRandom3 rnd2;
    double anglePhi1 = rnd1.Rndm()*TMath::Pi();
    double anglePhi2 = rnd2.Rndm()*TMath::Pi();

    DLM_Random RanGen(11);
    double RanVal1, RanVal2;

    TFile* F_EposDisto_p_pReso = new TFile(PipiIO::LegacyPRFileName);
    TNtuple* T_EposDisto_p_pReso = (TNtuple*)F_EposDisto_p_pReso->Get("InfoTuple_ClosePairs");
    if(!T_EposDisto_p_pReso){
        printf("ERROR: cannot find TNtuple InfoTuple_ClosePairs in ForMax_piReso_piReso_noOmega.root\n");
        delete F_EposDisto_p_pReso;
        return CleverMcLevyResoTM_pipi;
    }

    unsigned N_EposDisto_p_pReso = T_EposDisto_p_pReso->GetEntries();

    T_EposDisto_p_pReso->SetBranchAddress("k_D",&k_D);
    T_EposDisto_p_pReso->SetBranchAddress("P1",&fP1);
    T_EposDisto_p_pReso->SetBranchAddress("P2",&fP2);
    T_EposDisto_p_pReso->SetBranchAddress("M1",&fM1);
    T_EposDisto_p_pReso->SetBranchAddress("M2",&fM2);
    T_EposDisto_p_pReso->SetBranchAddress("Tau1",&Tau1);
    T_EposDisto_p_pReso->SetBranchAddress("Tau2",&Tau2);
    T_EposDisto_p_pReso->SetBranchAddress("AngleRcP1",&AngleRcP1);
    T_EposDisto_p_pReso->SetBranchAddress("AngleRcP2",&AngleRcP2);
    T_EposDisto_p_pReso->SetBranchAddress("AngleP1P2",&AngleP1P2);

    for(unsigned uEntry=0; uEntry<N_EposDisto_p_pReso; uEntry++){
        T_EposDisto_p_pReso->GetEntry(uEntry);
        Tau1 = 0;
        Tau2 = taureso;
        if(k_D>cutoff) continue;
        RanVal1 = (1+v2_cent*TMath::Cos(2*anglePhi1)) * RanGen.Exponential(fM2/(fP2*Tau2));
        CleverMcLevyResoTM_pipi->AddBGT_PR(RanVal1,-cos(AngleRcP2));
        CleverMcLevyResoTM_pipi->AddBGT_RP(RanVal1,cos(AngleRcP2));
    }
    delete F_EposDisto_p_pReso;

    TFile* F_EposDisto_pReso_pReso = new TFile(PipiIO::LegacyRRFileName);
    TNtuple* T_EposDisto_pReso_pReso = (TNtuple*)F_EposDisto_pReso_pReso->Get("InfoTuple_ClosePairs");
    if(!T_EposDisto_pReso_pReso){
        printf("ERROR: cannot find TNtuple InfoTuple_ClosePairs in EposDisto_pReso_pReso.root\n");
        delete F_EposDisto_pReso_pReso;
        return CleverMcLevyResoTM_pipi;
    }

    unsigned N_EposDisto_pReso_pReso = T_EposDisto_pReso_pReso->GetEntries();

    T_EposDisto_pReso_pReso->SetBranchAddress("k_D",&k_D);
    T_EposDisto_pReso_pReso->SetBranchAddress("P1",&fP1);
    T_EposDisto_pReso_pReso->SetBranchAddress("P2",&fP2);
    T_EposDisto_pReso_pReso->SetBranchAddress("M1",&fM1);
    T_EposDisto_pReso_pReso->SetBranchAddress("M2",&fM2);
    T_EposDisto_pReso_pReso->SetBranchAddress("Tau1",&Tau1);
    T_EposDisto_pReso_pReso->SetBranchAddress("Tau2",&Tau2);
    T_EposDisto_pReso_pReso->SetBranchAddress("AngleRcP1",&AngleRcP1);
    T_EposDisto_pReso_pReso->SetBranchAddress("AngleRcP2",&AngleRcP2);
    T_EposDisto_pReso_pReso->SetBranchAddress("AngleP1P2",&AngleP1P2);

    for(unsigned uEntry=0; uEntry<N_EposDisto_pReso_pReso; uEntry++){
        T_EposDisto_pReso_pReso->GetEntry(uEntry);

        if(k_D > cutoff) continue;
        if(fM1<=0.f || fM2<=0.f) continue;

        if(fP1<=0.0 || fP2<=0.0 || taureso<=0.0) continue;
        
        RanVal1 = (1+v2_cent*TMath::Cos(2*anglePhi1)) * RanGen.Exponential(fM1/(fP1*taureso));
        RanVal2 = (1+v2_cent*TMath::Cos(2*anglePhi2)) * RanGen.Exponential(fM2/(fP2*taureso));

        CleverMcLevyResoTM_pipi->AddBGT_RR(
            RanVal1, cos(AngleRcP1),
            RanVal2, cos(AngleRcP2),
            cos(AngleP1P2)
        );
    }
    delete F_EposDisto_pReso_pReso;

    if(FixAlpha && (Alpha==1 || Alpha==2)) CleverMcLevyResoTM_pipi->InitNumMcIter(1000000);
    else                                   CleverMcLevyResoTM_pipi->InitNumMcIter(100000);

    Kitty.SetAnaSource(CatsSourceForwarder, CleverMcLevyResoTM_pipi, 2);
    Kitty.SetAnaSource(0, SourceSize);
    Kitty.SetAnaSource(1, Alpha);

    Kitty.SetUseAnalyticSource(true);
    Kitty.SetAutoNormSource(false);

    return CleverMcLevyResoTM_pipi;
}

DLM_CleverMcDoubleGaussResoTM* SetDoubleGaussianResonanceSource(
    CATS& kitty,
    const double sourceCoreSize,
    const double sourceTailSize,
    const double sourceCoreFraction,
    const double shortLivedResonanceFraction,
    const double resonanceTau,
    const double cutoff)
{
  DLM_CleverMcDoubleGaussResoTM* source = new DLM_CleverMcDoubleGaussResoTM();
  source->InitScale(120, 2.0, 30.0);
  source->InitTailScale(120, 5.0, 160.0);
  source->InitCoreFraction(80, 0.0, 1.0);
  source->InitRad(500, 1e-4, 220.0);
  source->InitType(1);

  source->SetUpReso(0, shortLivedResonanceFraction);
  source->SetUpReso(1, shortLivedResonanceFraction);

  LoadOneDimensionalResonanceKernels(source,
                                     resonanceTau,
                                     cutoff);
  source->InitNumMcIter(1000000);

  kitty.SetAnaSource(CatsSourceForwarder, source, 3);
  kitty.SetAnaSource(0, sourceCoreSize);
  kitty.SetAnaSource(1, sourceTailSize);
  kitty.SetAnaSource(2, sourceCoreFraction);
  kitty.SetUseAnalyticSource(true);
  kitty.SetAutoNormSource(false);

  return source;
}

TTree* GetFirstExistingTree(TFile* file, const std::vector<TString>& treeNames)
{
    if(!file) return nullptr;
    for(const TString& name : treeNames){
        TTree* tree = dynamic_cast<TTree*>(file->Get(name));
        if(tree) return tree;
    }
    return nullptr;
}

TString GetFirstExistingBranch(TTree* tree, const std::vector<TString>& branchNames)
{
    if(!tree) return "";
    for(const TString& name : branchNames){
        if(tree->GetBranch(name)) return name;
    }
    return "";
}

bool SetVectorBranch(TTree* tree,
                     const std::vector<TString>& branchNames,
                     std::vector<double>*& values,
                     const char* label)
{
    const TString branchName = GetFirstExistingBranch(tree, branchNames);
    if(branchName==""){
        std::cout << "ERROR: cannot find branch for " << label
                  << " in tree " << (tree ? tree->GetName() : "null") << "\n";
        return false;
    }
    tree->SetBranchAddress(branchName, &values);
    return true;
}

void LoadPRFPrimaryResonancePairs(DLM_CleverMcLevyResoTM_3DTo1D* source,
                                  const TString& fileName,
                                  const double& tauReso)
{
    if(!source) return;

    TFile* file = TFile::Open(fileName, "READ");
    if(!file || file->IsZombie()){
        std::cout << "ERROR: cannot open PR/RP 3D EPOS file: " << fileName << "\n";
        if(file){ file->Close(); delete file; }
        return;
    }

    TTree* tree = GetFirstExistingTree(file, {"EPOS_PrimReson_PRF"});
    if(!tree){
        std::cout << "ERROR: cannot find EPOS_PrimReson_PRF in " << fileName << "\n";
        file->Close(); delete file;
        return;
    }

    std::vector<double>* pOut2 = nullptr;
    std::vector<double>* pSide2 = nullptr;
    std::vector<double>* pLong2 = nullptr;
    std::vector<double>* m2 = nullptr;
    std::vector<double>* oS_dot_oD = nullptr;
    std::vector<double>* sS_dot_oD = nullptr;
    std::vector<double>* oS_dot_sD = nullptr;
    std::vector<double>* sS_dot_sD = nullptr;
    std::vector<double>* longScale = nullptr;

    const bool ok =
        SetVectorBranch(tree, {"motherInDaughterPRF_pOut2"}, pOut2, "mother pOut2") &&
        SetVectorBranch(tree, {"motherInDaughterPRF_pSide2"}, pSide2, "mother pSide2") &&
        SetVectorBranch(tree, {"motherInDaughterPRF_pLong2"}, pLong2, "mother pLong2") &&
        SetVectorBranch(tree, {"m2"}, m2, "m2") &&
        SetVectorBranch(tree, {"oS_dot_oD"}, oS_dot_oD, "oS_dot_oD") &&
        SetVectorBranch(tree, {"sS_dot_oD"}, sS_dot_oD, "sS_dot_oD") &&
        SetVectorBranch(tree, {"oS_dot_sD"}, oS_dot_sD, "oS_dot_sD") &&
        SetVectorBranch(tree, {"sS_dot_sD"}, sS_dot_sD, "sS_dot_sD") &&
        SetVectorBranch(tree, {"longScale"}, longScale, "longScale");

    if(!ok){
        file->Close(); delete file;
        return;
    }

    unsigned nAdded = 0;
    for(Long64_t iEntry=0; iEntry<tree->GetEntries(); ++iEntry){
        tree->GetEntry(iEntry);
        if(!pOut2 || !pSide2 || !pLong2 || !m2 ||
           !oS_dot_oD || !sS_dot_oD || !oS_dot_sD || !sS_dot_sD || !longScale) continue;

        const size_t nPairs = std::min({
            pOut2->size(), pSide2->size(), pLong2->size(), m2->size(),
            oS_dot_oD->size(), sS_dot_oD->size(), oS_dot_sD->size(),
            sS_dot_sD->size(), longScale->size()
        });

        for(size_t iPair=0; iPair<nPairs; ++iPair){
            const double mass = m2->at(iPair);
            if(mass<=0.) continue;

            const TVector3 motherMomentum(
                pOut2->at(iPair),
                pSide2->at(iPair),
                pLong2->at(iPair)
            );

            source->AddBGT3D_PR(
                motherMomentum, mass, tauReso,
                oS_dot_oD->at(iPair), sS_dot_oD->at(iPair),
                oS_dot_sD->at(iPair), sS_dot_sD->at(iPair),
                longScale->at(iPair)
            );
            source->AddBGT3D_RP(
                motherMomentum, mass, tauReso,
                oS_dot_oD->at(iPair), sS_dot_oD->at(iPair),
                oS_dot_sD->at(iPair), sS_dot_sD->at(iPair),
                longScale->at(iPair)
            );
            ++nAdded;
        }
    }

    std::cout << "Loaded " << nAdded << " PR and RP 3D PRF kernels from "
              << fileName << " (" << tree->GetName() << ")\n";
    file->Close();
    delete file;
}

void LoadPRFResonanceResonancePairs(DLM_CleverMcLevyResoTM_3DTo1D* source,
                                    const TString& fileName,
                                    const double& tauReso)
{
    if(!source) return;

    TFile* file = TFile::Open(fileName, "READ");
    if(!file || file->IsZombie()){
        std::cout << "ERROR: cannot open RR 3D EPOS file: " << fileName << "\n";
        if(file){ file->Close(); delete file; }
        return;
    }

    TTree* tree = GetFirstExistingTree(file, {"EPOS_Reson_PRF"});
    if(!tree){
        std::cout << "ERROR: cannot find EPOS_Reson_PRF in " << fileName << "\n";
        file->Close(); delete file;
        return;
    }

    std::vector<double>* pOut1 = nullptr;
    std::vector<double>* pSide1 = nullptr;
    std::vector<double>* pLong1 = nullptr;
    std::vector<double>* m1 = nullptr;
    std::vector<double>* pOut2 = nullptr;
    std::vector<double>* pSide2 = nullptr;
    std::vector<double>* pLong2 = nullptr;
    std::vector<double>* m2 = nullptr;
    std::vector<double>* oS_dot_oD = nullptr;
    std::vector<double>* sS_dot_oD = nullptr;
    std::vector<double>* oS_dot_sD = nullptr;
    std::vector<double>* sS_dot_sD = nullptr;
    std::vector<double>* longScale = nullptr;

    const bool ok =
        SetVectorBranch(tree, {"motherInDaughterPRF_pOut1"}, pOut1, "mother pOut1") &&
        SetVectorBranch(tree, {"motherInDaughterPRF_pSide1"}, pSide1, "mother pSide1") &&
        SetVectorBranch(tree, {"motherInDaughterPRF_pLong1"}, pLong1, "mother pLong1") &&
        SetVectorBranch(tree, {"m1"}, m1, "m1") &&
        SetVectorBranch(tree, {"motherInDaughterPRF_pOut2"}, pOut2, "mother pOut2") &&
        SetVectorBranch(tree, {"motherInDaughterPRF_pSide2"}, pSide2, "mother pSide2") &&
        SetVectorBranch(tree, {"motherInDaughterPRF_pLong2"}, pLong2, "mother pLong2") &&
        SetVectorBranch(tree, {"m2"}, m2, "m2") &&
        SetVectorBranch(tree, {"oS_dot_oD"}, oS_dot_oD, "oS_dot_oD") &&
        SetVectorBranch(tree, {"sS_dot_oD"}, sS_dot_oD, "sS_dot_oD") &&
        SetVectorBranch(tree, {"oS_dot_sD"}, oS_dot_sD, "oS_dot_sD") &&
        SetVectorBranch(tree, {"sS_dot_sD"}, sS_dot_sD, "sS_dot_sD") &&
        SetVectorBranch(tree, {"longScale"}, longScale, "longScale");

    if(!ok){
        file->Close(); delete file;
        return;
    }

    unsigned nAdded = 0;
    for(Long64_t iEntry=0; iEntry<tree->GetEntries(); ++iEntry){
        tree->GetEntry(iEntry);
        if(!pOut1 || !pSide1 || !pLong1 || !m1 ||
           !pOut2 || !pSide2 || !pLong2 || !m2 ||
           !oS_dot_oD || !sS_dot_oD || !oS_dot_sD || !sS_dot_sD || !longScale) continue;

        const size_t nPairs = std::min({
            pOut1->size(), pSide1->size(), pLong1->size(), m1->size(),
            pOut2->size(), pSide2->size(), pLong2->size(), m2->size(),
            oS_dot_oD->size(), sS_dot_oD->size(), oS_dot_sD->size(),
            sS_dot_sD->size(), longScale->size()
        });

        for(size_t iPair=0; iPair<nPairs; ++iPair){
            const double mass1 = m1->at(iPair);
            const double mass2 = m2->at(iPair);
            if(mass1<=0. || mass2<=0.) continue;

            const TVector3 motherMomentum1(
                pOut1->at(iPair),
                pSide1->at(iPair),
                pLong1->at(iPair)
            );
            const TVector3 motherMomentum2(
                pOut2->at(iPair),
                pSide2->at(iPair),
                pLong2->at(iPair)
            );

            source->AddBGT3D_RR(
                motherMomentum1, mass1, tauReso,
                motherMomentum2, mass2, tauReso,
                oS_dot_oD->at(iPair), sS_dot_oD->at(iPair),
                oS_dot_sD->at(iPair), sS_dot_sD->at(iPair),
                longScale->at(iPair)
            );
            ++nAdded;
        }
    }

    std::cout << "Loaded " << nAdded << " RR 3D PRF kernels from "
              << fileName << " (" << tree->GetName() << ")\n";
    file->Close();
    delete file;
}

DLM_CleverMcLevyResoTM_3DTo1D* CATS_ResoSource_pipi_3DTo1D_UsingPionAngles(
    CATS& Kitty,
    const double& Rout,
    const double& Rside,
    const double& Rlong,
    const double& Alpha,
    const bool& FixAlpha,
    const double& fracreso,
    const double& taureso,
    const double& massreso,
    const double& cutoff,
    const bool& usePRF3DShift=false
){
    (void)massreso;

    DLM_CleverMcLevyResoTM_3DTo1D* CleverMcLevyResoTM_pipi =
        new DLM_CleverMcLevyResoTM_3DTo1D();

    if(FixAlpha) CleverMcLevyResoTM_pipi->InitStability(1, Alpha-1e-6, Alpha+1e-6);
    else         CleverMcLevyResoTM_pipi->InitStability(21, 1, 2);

    CleverMcLevyResoTM_pipi->InitScale(100, 3.0, 20.0);
    CleverMcLevyResoTM_pipi->InitRad(400, 1e-4, 128);
    CleverMcLevyResoTM_pipi->InitType(2);
    CleverMcLevyResoTM_pipi->SetShiftMode(
        usePRF3DShift ?
        DLM_CleverMcLevyResoTM_3DTo1D::rsmVector3D :
        DLM_CleverMcLevyResoTM_3DTo1D::rsmLegacy1D
    );

    CleverMcLevyResoTM_pipi->SetUpReso(0, fracreso);
    CleverMcLevyResoTM_pipi->SetUpReso(1, fracreso);

    if(usePRF3DShift){
        LoadPRFPrimaryResonancePairs(
            CleverMcLevyResoTM_pipi,
            PipiIO::PRF3DPRFileName,
            taureso
        );
        LoadPRFResonanceResonancePairs(
            CleverMcLevyResoTM_pipi,
            PipiIO::PRF3DRRFileName,
            taureso
        );

        if(FixAlpha && (Alpha==1 || Alpha==2)) CleverMcLevyResoTM_pipi->InitNumMcIter(1000000);
        else                                   CleverMcLevyResoTM_pipi->InitNumMcIter(100000);

        Kitty.SetAnaSource(CatsSourceForwarder, CleverMcLevyResoTM_pipi, 4);
        Kitty.SetAnaSource(0, Rout);
        Kitty.SetAnaSource(1, Rside);
        Kitty.SetAnaSource(2, Rlong);
        Kitty.SetAnaSource(3, Alpha);

        Kitty.SetUseAnalyticSource(true);
        Kitty.SetAutoNormSource(false);

        return CleverMcLevyResoTM_pipi;
    }

    Float_t k_D, fP1, fP2, fM1, fM2, Tau1, Tau2;
    Float_t AngleRcP1, AngleRcP2, AngleP1P2;

    DLM_Random RanGen(11);
    double RanVal1, RanVal2;

    TFile* F_EposDisto_p_pReso = new TFile(PipiIO::LegacyPRFileName);
    TNtuple* T_EposDisto_p_pReso = (TNtuple*)F_EposDisto_p_pReso->Get("InfoTuple_ClosePairs");
    if(!T_EposDisto_p_pReso){
        printf("ERROR: cannot find TNtuple InfoTuple_ClosePairs in ForMax_pi_piReso_noOmega.root\n");
        delete F_EposDisto_p_pReso;
        return CleverMcLevyResoTM_pipi;
    }

    unsigned N_EposDisto_p_pReso = T_EposDisto_p_pReso->GetEntries();

    T_EposDisto_p_pReso->SetBranchAddress("k_D",&k_D);
    T_EposDisto_p_pReso->SetBranchAddress("P1",&fP1);
    T_EposDisto_p_pReso->SetBranchAddress("P2",&fP2);
    T_EposDisto_p_pReso->SetBranchAddress("M1",&fM1);
    T_EposDisto_p_pReso->SetBranchAddress("M2",&fM2);
    T_EposDisto_p_pReso->SetBranchAddress("Tau1",&Tau1);
    T_EposDisto_p_pReso->SetBranchAddress("Tau2",&Tau2);
    T_EposDisto_p_pReso->SetBranchAddress("AngleRcP1",&AngleRcP1);
    T_EposDisto_p_pReso->SetBranchAddress("AngleRcP2",&AngleRcP2);
    T_EposDisto_p_pReso->SetBranchAddress("AngleP1P2",&AngleP1P2);

    for(unsigned uEntry=0; uEntry<N_EposDisto_p_pReso; uEntry++){
        T_EposDisto_p_pReso->GetEntry(uEntry);
        Tau1 = 0;
        Tau2 = taureso;
        if(k_D>cutoff) continue;
        RanVal1 = RanGen.Exponential(fM2/(fP2*Tau2));
        CleverMcLevyResoTM_pipi->AddBGT_PR(RanVal1,-cos(AngleRcP2));
        CleverMcLevyResoTM_pipi->AddBGT_RP(RanVal1,cos(AngleRcP2));
    }
    delete F_EposDisto_p_pReso;

    TFile* F_EposDisto_pReso_pReso = new TFile(PipiIO::LegacyRRFileName);
    TNtuple* T_EposDisto_pReso_pReso = (TNtuple*)F_EposDisto_pReso_pReso->Get("InfoTuple_ClosePairs");
    if(!T_EposDisto_pReso_pReso){
        printf("ERROR: cannot find TNtuple InfoTuple_ClosePairs in ForMax_piReso_piReso_noOmega.root\n");
        delete F_EposDisto_pReso_pReso;
        return CleverMcLevyResoTM_pipi;
    }

    unsigned N_EposDisto_pReso_pReso = T_EposDisto_pReso_pReso->GetEntries();

    T_EposDisto_pReso_pReso->SetBranchAddress("k_D",&k_D);
    T_EposDisto_pReso_pReso->SetBranchAddress("P1",&fP1);
    T_EposDisto_pReso_pReso->SetBranchAddress("P2",&fP2);
    T_EposDisto_pReso_pReso->SetBranchAddress("M1",&fM1);
    T_EposDisto_pReso_pReso->SetBranchAddress("M2",&fM2);
    T_EposDisto_pReso_pReso->SetBranchAddress("Tau1",&Tau1);
    T_EposDisto_pReso_pReso->SetBranchAddress("Tau2",&Tau2);
    T_EposDisto_pReso_pReso->SetBranchAddress("AngleRcP1",&AngleRcP1);
    T_EposDisto_pReso_pReso->SetBranchAddress("AngleRcP2",&AngleRcP2);
    T_EposDisto_pReso_pReso->SetBranchAddress("AngleP1P2",&AngleP1P2);

    for(unsigned uEntry=0; uEntry<N_EposDisto_pReso_pReso; uEntry++){
        T_EposDisto_pReso_pReso->GetEntry(uEntry);
        if(k_D > cutoff) continue;
        if(fM1<=0.f || fM2<=0.f) continue;
        if(fP1<=0.0 || fP2<=0.0 || taureso<=0.0) continue;

        RanVal1 = RanGen.Exponential(fM1/(fP1*taureso));
        RanVal2 = RanGen.Exponential(fM2/(fP2*taureso));

        CleverMcLevyResoTM_pipi->AddBGT_RR(
            RanVal1, cos(AngleRcP1),
            RanVal2, cos(AngleRcP2),
            cos(AngleP1P2)
        );
    }
    delete F_EposDisto_pReso_pReso;

    if(FixAlpha && (Alpha==1 || Alpha==2)) CleverMcLevyResoTM_pipi->InitNumMcIter(1000000);
    else                                   CleverMcLevyResoTM_pipi->InitNumMcIter(100000);

    Kitty.SetAnaSource(CatsSourceForwarder, CleverMcLevyResoTM_pipi, 4);
    Kitty.SetAnaSource(0, Rout);
    Kitty.SetAnaSource(1, Rside);
    Kitty.SetAnaSource(2, Rlong);
    Kitty.SetAnaSource(3, Alpha);

    Kitty.SetUseAnalyticSource(true);
    Kitty.SetAutoNormSource(false);

    return CleverMcLevyResoTM_pipi;
}

// =========================================================
// Helper: basic pi-pi settings (QS + optional Coulomb, no strong)
// =========================================================
void CATS_pipi_Basic(CATS& Kitty,
                     const unsigned NumMomBins,
                     const double kMin,
                     const double kMax,
                     const bool Identical=true,
                     const bool WithCoulomb=true){

  Kitty.SetMomBins(NumMomBins, kMin, kMax);

  Kitty.SetNumChannels(1);
  Kitty.SetNumPW(0,0);
  Kitty.SetSpin(0,0);
  Kitty.SetChannelWeight(0, 1.0);

  // like-sign pi pi: Q1Q2 = +1 (only meaningful if Coulomb on)
  if(WithCoulomb) Kitty.SetQ1Q2(1);
  else            Kitty.SetQ1Q2(0);

  if(Identical) Kitty.SetPdgId(211, 211);
  else          Kitty.SetPdgId(-211, 211);

  if(WithCoulomb)      Kitty.SetRedMass(0.5*Mass_Pic);
  else if(Identical)   Kitty.SetRedMass(0.5*Mass_Pi0);
  else                 Kitty.SetRedMass((Mass_Pi0*Mass_Pic)/(Mass_Pi0+Mass_Pic));
}

DLM_Histo<double>* MakeDlmFeedDownFractionFromLambdaTH1(const TH1* hLambda,
                                                        const double xScale=1000.0)
{
  if(!hLambda) return nullptr;

  const unsigned nBins = hLambda->GetNbinsX();
  std::vector<double> edges(nBins+1);
  for(unsigned iBin=0; iBin<=nBins; ++iBin){
    edges[iBin] = hLambda->GetXaxis()->GetBinLowEdge(iBin+1) * xScale;
  }

  DLM_Histo<double>* hFeedDown = new DLM_Histo<double>();
  hFeedDown->SetUp(1);
  hFeedDown->SetUp(0, nBins, edges.data());
  hFeedDown->Initialize(true);

  unsigned whichBin[1];
  for(unsigned iBin=0; iBin<nBins; ++iBin){
    const double lambdaMain = hLambda->GetBinContent(iBin+1) * (1-fracWeakDecay) * (1-fracWeakDecay);
    const double lambdaErr = hLambda->GetBinError(iBin+1) * (1-fracWeakDecay) * (1-fracWeakDecay);

    if(!std::isfinite(lambdaMain) || lambdaMain<0.0 || lambdaMain>1.0){
      std::cout << "ERROR: lambda histogram has invalid value in bin "
                << iBin+1 << ": " << lambdaMain << "\n";
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
    return nullptr;
  }

  TH1* hLambda = dynamic_cast<TH1*>(file->Get(histName));
  if(!hLambda){
    std::cout << "ERROR: cannot find lambda histogram: " << histName << "\n";
    file->Close();
    return nullptr;
  }

  DLM_Histo<double>* hFeedDown = MakeDlmFeedDownFractionFromLambdaTH1(hLambda, xScale);
  file->Close();
  return hFeedDown;
}

// =========================================================
// Init theory objects (call once before Fit)
// =========================================================
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
  // cleanup old objects if re-initializing in same session
  if(gCkPiPi){ delete gCkPiPi; gCkPiPi = nullptr; }
  if(gKittyPiPi){ delete gKittyPiPi; gKittyPiPi = nullptr; }
  if(gCRSpipi){ delete gCRSpipi; gCRSpipi=nullptr; }
  if(gCRSpipi3DTo1D){ delete gCRSpipi3DTo1D; gCRSpipi3DTo1D=nullptr; }
  if(gCRSpipiDoubleGaus){ delete gCRSpipiDoubleGaus; gCRSpipiDoubleGaus=nullptr; }
  if(gFemto_pipi_fitter){ delete gFemto_pipi_fitter; gFemto_pipi_fitter=nullptr; }
  gUseDoubleGaussianSource = false;

  gKittyPiPi = new CATS();

  if(SourceType=="Gauss"){
    CATS_GaussSource(*gKittyPiPi, SourceSize);
  }
  else if(SourceType=="CoreReso3DTo1D" || SourceType=="CoreReso3DTo1D_PRF"){
    const bool usePRF3DShift = (SourceType=="CoreReso3DTo1D_PRF");
    gCRSpipi3DTo1D = CATS_ResoSource_pipi_3DTo1D_UsingPionAngles(
        *gKittyPiPi,
        SourceSize, SourceSize, SourceSize,
        SourceAlpha, true,
        fracShortLivedResonance,
        1.5,   // taureso [fm]
        1124,  // massreso [MeV]
        kMax,  // cutoff [MeV]
        usePRF3DShift
    );
  }else if(SourceType=="CoreResoDoubleGauss"){
    gCRSpipiDoubleGaus = SetDoubleGaussianResonanceSource(
        *gKittyPiPi, 
        SourceSize,
        doubleGaussTailWidth,
        doubleGaussCoreFraction,
        fracShortLivedResonance,
        1.5,   // taureso [fm]
        kMax   // cutoff [MeV]    
    );
    gUseDoubleGaussianSource = true;
  }else{
    gCRSpipi = CATS_ResoSource_pipi_UsingPionAngles(
        *gKittyPiPi, SourceSize,
        SourceAlpha, true,
        fracShortLivedResonance,
        1.5,   // taureso [fm]
        1124,  // massreso [MeV]
        kMax   // cutoff [MeV]
    );
  }

  CATS_pipi_Basic(*gKittyPiPi, NumMomBins, kMin, kMax, Identical, WithCoulomb);

  // compute once to allocate internals
  gKittyPiPi->SetNotifications(CATS::nWarning);
  gKittyPiPi->SetNotifications(CATS::nError);
  gKittyPiPi->KillTheCat();

  // DLM_Ck wrapper around CATS
  gCkPiPi = new DLM_Ck(gKittyPiPi->GetNumSourcePars(), 0, *gKittyPiPi);
  gCkPiPi->Update(true);

  TH2F* hMomentumResolution_pipi =
      GetSmearMatrix(PipiIO::SmearingFileName, PipiIO::SmearingHistoName);
  if(!hMomentumResolution_pipi){
    std::cout << "ERROR: failed to get Smear Matrix\n";
    return;
  }

  DLM_Histo<double>* hLambdaFeedDown =
      LoadDlmFeedDownFraction(lambdaFileName, lambdaHistName, 1000.0);
  if(!hLambdaFeedDown){
    std::cout << "ERROR: failed to build k*-dependent lambda for AddContribution\n";
    return;
  }

  DLM_CkDecomposition* CkDec_pipi = new DLM_CkDecomposition("pipi",1,*gCkPiPi,hMomentumResolution_pipi);
  CkDec_pipi->AddContribution(0,*hLambdaFeedDown,DLM_CkDecomposition::cFeedDown);

  delete hLambdaFeedDown;

  CkDec_pipi->Update(true,true);
  gFemto_pipi_fitter = CkDec_pipi;
}

double femtoFunc(double* x, double* par){
  double& MOM = x[0];
  double& source_size = par[0];
  double& source_alpha = par[1];

  static double lastR = 1e99;
  static double lastAlpha = 1e99;
  static bool first = true;

  // Update the CATS source only when the source parameters actually change.
  if(first || fabs(source_size-lastR) > 1e-6 || fabs(source_alpha-lastAlpha) > 1e-6){
    lastR = par[0];
    lastAlpha = par[1];
    first = false;

    gFemto_pipi_fitter->GetCk()->SetSourcePar(0,source_size);
    if(gUseDoubleGaussianSource){
      gFemto_pipi_fitter->GetCk()->SetSourcePar(1,doubleGaussTailWidth);
      gFemto_pipi_fitter->GetCk()->SetSourcePar(2,doubleGaussCoreFraction);
    }else{
      gFemto_pipi_fitter->GetCk()->SetSourcePar(1,source_alpha);
    }
    gFemto_pipi_fitter->Update(true);
  }

  return gFemto_pipi_fitter->EvalCk(MOM);
}

double femtoFunc3DTo1D(double* x, double* par){
  double& MOM = x[0];
  double& source_rout = par[0];
  double& source_rside = par[1];
  double& source_rlong = par[2];
  double& source_alpha = par[3];

  static double lastRout = 1e99;
  static double lastRside = 1e99;
  static double lastRlong = 1e99;
  static double lastAlpha = 1e99;
  static bool first = true;

  if(first ||
     fabs(source_rout-lastRout) > 1e-6 ||
     fabs(source_rside-lastRside) > 1e-6 ||
     fabs(source_rlong-lastRlong) > 1e-6 ||
     fabs(source_alpha-lastAlpha) > 1e-6){
    lastRout = par[0];
    lastRside = par[1];
    lastRlong = par[2];
    lastAlpha = par[3];
    first = false;

    gFemto_pipi_fitter->GetCk()->SetSourcePar(0,source_rout);
    gFemto_pipi_fitter->GetCk()->SetSourcePar(1,source_rside);
    gFemto_pipi_fitter->GetCk()->SetSourcePar(2,source_rlong);
    gFemto_pipi_fitter->GetCk()->SetSourcePar(3,source_alpha);
    gFemto_pipi_fitter->Update(true);
  }

  return gFemto_pipi_fitter->EvalCk(MOM);
}

// =========================================================
// TF1 callback: C(k) * baseline
//
// par layout:
// 0: R, Levy scale/source size [fm]
// 1: alpha, Levy stability; alpha=2 is Gaussian, alpha<2 gives a longer tail
// 2: baseline p0
// 3: baseline p1
// =========================================================
double MickeyFitter_pipi_Pol3(double* x, double* par)
{
  const double MOM = x[0];

  if(!gKittyPiPi || !gCkPiPi || !gFemto_pipi_fitter){
    // If user forgot to call Init_pipi_Theory, avoid segfault and return 1
    return 1.0;
  }

  const double Ck_core = femtoFunc(x,&par[0]);

  const double Ck_tot = Ck_core;

  double &p0 = par[2];
  double &p1 = par[3];
  const double BaseLine = p0 * (1. + p1*MOM );

  return  Ck_tot * BaseLine;
}

// par layout:
// 0: Rout, 1: Rside, 2: Rlong, 3: alpha, 4: baseline p0, 5: baseline p1
double MickeyFitter_pipi_3DTo1D_Pol1(double* x, double* par)
{
  const double MOM = x[0];

  if(!gKittyPiPi || !gCkPiPi || !gFemto_pipi_fitter){
    return 1.0;
  }

  const double Ck_tot = femtoFunc3DTo1D(x,&par[0]);
  double &p0 = par[4];
  double &p1 = par[5];
  const double BaseLine = p0 * (1. + p1*MOM);

  return Ck_tot * BaseLine;
}

double baseline3(double* x, double* p){
  const double xx = x[0];
  return p[0] * (1.0 + p[1]*xx + p[2]*xx*xx + p[3]*xx*xx*xx);
}

TH1D* CloneWithScaledXAxis(const TH1* hIn, const char* newName, const double xScale)
{
  if(!hIn) return nullptr;

  const int n = hIn->GetNbinsX();

  std::vector<double> edges(n+1);
  for(int i=0;i<=n;i++){
    edges[i] = hIn->GetXaxis()->GetBinLowEdge(i+1) * xScale;
  }

  TH1D* hOut = new TH1D(newName, hIn->GetTitle(), n, edges.data());
  hOut->Sumw2();

  for(int i=1;i<=n;i++){
    hOut->SetBinContent(i, hIn->GetBinContent(i));
    hOut->SetBinError(i,   hIn->GetBinError(i));
  }

  hOut->GetXaxis()->SetTitle("k (MeV/c)");
  return hOut;
}

TH1D* FillHistoFromTF1(const char* hname,
                       TF1* func,
                       int nbins,
                       double xmin,
                       double xmax)
{
    TH1D* h = new TH1D(hname, "", nbins, xmin, xmax);

    for (int i = 1; i <= nbins; i++) {
        double x = h->GetBinCenter(i);
        double y = func->Eval(x);
        h->SetBinContent(i, y);
    }

    return h;
}

TH1D* RebinCfWeighted(const TH1* hIn,
                      const char* newName,
                      const std::vector<double>& edges)
{
  if(!hIn || edges.size()<2) return nullptr;

  TH1D* hOut = new TH1D(newName, hIn->GetTitle(), (int)edges.size()-1, edges.data());
  hOut->Sumw2();
  hOut->GetXaxis()->SetTitle(hIn->GetXaxis()->GetTitle());
  hOut->GetYaxis()->SetTitle(hIn->GetYaxis()->GetTitle());

  for(int iOut=1; iOut<=hOut->GetNbinsX(); ++iOut){
    const double xLow = hOut->GetXaxis()->GetBinLowEdge(iOut);
    const double xUp = hOut->GetXaxis()->GetBinUpEdge(iOut);
    double sumW = 0.0;
    double sumWY = 0.0;
    int nUsed = 0;

    for(int iIn=1; iIn<=hIn->GetNbinsX(); ++iIn){
      const double x = hIn->GetXaxis()->GetBinCenter(iIn);
      if(x<xLow || x>=xUp) continue;

      const double y = hIn->GetBinContent(iIn);
      const double ey = hIn->GetBinError(iIn);
      if(!std::isfinite(y)) continue;

      const double w = (ey>0.0 && std::isfinite(ey)) ? 1.0/(ey*ey) : 1.0;
      sumW += w;
      sumWY += w*y;
      ++nUsed;
    }

    if(sumW>0.0 && nUsed>0){
      hOut->SetBinContent(iOut, sumWY/sumW);
      hOut->SetBinError(iOut, std::sqrt(1.0/sumW));
    }
    else{
      hOut->SetBinContent(iOut, 0.0);
      hOut->SetBinError(iOut, 0.0);
    }
  }

  return hOut;
}

std::vector<double> UniformEdges(double xMin, double xMax, double binWidth)
{
  std::vector<double> edges;
  if(binWidth<=0.0 || xMax<=xMin) return edges;

  for(double x=xMin; x<xMax-1e-9; x+=binWidth){
    edges.push_back(x);
  }
  if(edges.empty() || fabs(edges.back()-xMax)>1e-9) edges.push_back(xMax);
  return edges;
}

std::vector<double> VariablePbPbFitEdges()
{
  std::vector<double> edges;
  for(double x=0.0; x<24.0-1e-9; x+=4.0) edges.push_back(x);
  for(double x=24.0; x<80.0-1e-9; x+=8.0) edges.push_back(x);
  for(double x=80.0; x<200.0-1e-9; x+=16.0) edges.push_back(x);
  edges.push_back(200.0);
  return edges;
}

TF1* FitFixedAlphaOnHistogram(TH1* hData,
                              const char* fitName,
                              const double fixedAlpha,
                              const double sourceSizeInit,
                              const double fitMin,
                              const double fitMax,
                              const double baselineFitMin,
                              const double baselineFitMax)
{
  if(!hData) return nullptr;

  TF1* baseline_prefit = new TF1(Form("%s_baseline_prefit", fitName),
                                 "[0]*(1.0+[1]*x)", baselineFitMin, baselineFitMax);
  baseline_prefit->SetParameter(0, 1.0);
  baseline_prefit->SetParameter(1, 0.0);
  hData->Fit(baseline_prefit, "RQ0", "", baselineFitMin, baselineFitMax);

  TF1* fit = new TF1(fitName, MickeyFitter_pipi_Pol3, fitMin, fitMax, 4);
  fit->SetNpx(PipiIO::FitNpx);

  fit->SetParameter(0, sourceSizeInit);
  fit->SetParLimits(0, 8.0, 18.0);
  fit->SetParameter(1, fixedAlpha);
  fit->FixParameter(1, fixedAlpha);
  fit->SetParameter(2, baseline_prefit->GetParameter(0));
  fit->FixParameter(2, baseline_prefit->GetParameter(0));
  fit->SetParameter(3, baseline_prefit->GetParameter(1));
  fit->FixParameter(3, baseline_prefit->GetParameter(1));

  hData->Fit(fit, "RSQ0", "", fitMin, fitMax);
  delete baseline_prefit;
  return fit;
}

void ScanAlphaBinning(bool Identical,
                      bool WithCoulomb,
                      double alphaMin,
                      double alphaMax,
                      double alphaStep,
                      double fitMin,
                      double fitMax,
                      double baselineFitMin,
                      double baselineFitMax,
                      unsigned catsMomBins)
{
  (void)fitMin;
  (void)fitMax;
  (void)baselineFitMin;
  (void)baselineFitMax;
  (void)catsMomBins;

  const TString outputFileName = PipiIO::OutputScanAlphaBinning;

  TFile* inputFile = TFile::Open(PipiIO::DataFileName, "READ");
  if(!inputFile || inputFile->IsZombie()){
    std::cout << "ERROR: cannot open input ROOT file: " << PipiIO::DataFileName << "\n";
    return;
  }

  TH1* hDataRaw = dynamic_cast<TH1*>(inputFile->Get(PipiIO::DataHistoName));
  if(!hDataRaw){
    std::cout << "ERROR: cannot find histogram: " << PipiIO::DataHistoName << "\n";
    inputFile->Close();
    return;
  }

  std::vector<std::pair<TString, std::vector<double>>> binnings;
  binnings.push_back({"bin4MeV", UniformEdges(0.0, 200.0, 4.0)});
  binnings.push_back({"bin8MeV", UniformEdges(0.0, 200.0, 8.0)});
  binnings.push_back({"bin12MeV", UniformEdges(0.0, 204.0, 12.0)});
  binnings.push_back({"binVariable_4_8_16MeV", VariablePbPbFitEdges()});

  TFile* outputFile = new TFile(outputFileName, "RECREATE");
  TNtuple summary("alpha_binning_scan",
                  "alpha_binning_scan",
                  "binning:alpha:R:chi2:ndf:chi2ndf:status:p0:p1");

  const double sourceSizeInit = 12.3;
  unsigned iBinning = 0;

  for(auto& binning : binnings){
    const TString binningName = binning.first;
    TH1D* hRebinned = RebinCfWeighted(hDataRaw, Form("data_%s", binningName.Data()), binning.second);
    if(!hRebinned) continue;

    outputFile->cd();
    hRebinned->Write();

    for(double alpha=alphaMin; alpha<=alphaMax+1e-9; alpha+=alphaStep){
      std::cout << "ScanAlphaBinning: " << binningName
                << ", alpha=" << alpha << "\n";

      Init_pipi_Theory("CoreReso",
                       Identical,
                       WithCoulomb,
                       0,
                       sourceSizeInit,
                       alpha,
                       PipiIO::CatsNumMomBins,
                       PipiIO::CatsKMin,
                       PipiIO::CatsKMax,
                       PipiIO::LambdaFileName,
                       PipiIO::LambdaHistoName);

      TF1* fit = FitFixedAlphaOnHistogram(hRebinned,
                                          Form("fit_%s_alpha%.2f", binningName.Data(), alpha),
                                          alpha,
                                          sourceSizeInit,
                                          PipiIO::FitMin,
                                          PipiIO::FitMax,
                                          PipiIO::BaselineFitMin,
                                          PipiIO::BaselineFitMax);
      if(!fit) continue;

      const double chi2 = fit->GetChisquare();
      const int ndf = fit->GetNDF();
      const double chi2ndf = ndf>0 ? chi2/double(ndf) : 0.0;
      summary.Fill(iBinning,
                   alpha,
                   fit->GetParameter(0),
                   chi2,
                   ndf,
                   chi2ndf,
                   0,
                   fit->GetParameter(2),
                   fit->GetParameter(3));

      outputFile->cd();
      fit->Write("", TObject::kOverwrite);
      summary.Write("", TObject::kOverwrite);
      outputFile->Save();
      delete fit;
    }

    ++iBinning;
  }

  outputFile->cd();
  summary.Write();
  outputFile->Close();
  inputFile->Close();
  delete outputFile;
  delete inputFile;

  std::cout << "Saved alpha/binning scan to " << outputFileName << "\n";
}

void ScanAlpha8MeV(bool Identical,
                   bool WithCoulomb,
                   double alphaMin,
                   double alphaMax,
                   double alphaStep,
                   double fitMin,
                   double fitMax,
                   double baselineFitMin,
                   double baselineFitMax,
                   unsigned catsMomBins)
{
  (void)fitMin;
  (void)fitMax;
  (void)baselineFitMin;
  (void)baselineFitMax;
  (void)catsMomBins;

  const double savedWeakDecay = fracWeakDecay;

  const TString outputFileName = PipiIO::OutputScanAlpha8MeV;

  TFile* inputFile = TFile::Open(PipiIO::DataFileName, "READ");
  if(!inputFile || inputFile->IsZombie()){
    std::cout << "ERROR: cannot open input ROOT file: " << PipiIO::DataFileName << "\n";
    fracWeakDecay = savedWeakDecay;
    return;
  }

  TH1* hDataRaw = dynamic_cast<TH1*>(inputFile->Get(PipiIO::DataHistoName));
  if(!hDataRaw){
    std::cout << "ERROR: cannot find histogram: " << PipiIO::DataHistoName << "\n";
    inputFile->Close();
    fracWeakDecay = savedWeakDecay;
    return;
  }

  TH1D* hRebinned = RebinCfWeighted(hDataRaw, "data_bin8MeV", UniformEdges(0.0, 200.0, 4.0));
  if(!hRebinned){
    inputFile->Close();
    fracWeakDecay = savedWeakDecay;
    return;
  }

  TFile* outputFile = new TFile(outputFileName, "RECREATE");
  TNtuple summary("alpha8MeV_scan",
                  "alpha8MeV_scan",
                  "alpha:R:chi2:ndf:chi2ndf:p0:p1");

  outputFile->cd();
  hRebinned->Write();

  const double sourceSizeInit = 12.3;
  for(double alpha=alphaMin; alpha<=alphaMax+1e-9; alpha+=alphaStep){
    std::cout << "ScanAlpha8MeV: alpha=" << alpha << "\n";

    Init_pipi_Theory("CoreReso",
                     Identical,
                     WithCoulomb,
                     0,
                     sourceSizeInit,
                     alpha,
                     PipiIO::CatsNumMomBins,
                     PipiIO::CatsKMin,
                     PipiIO::CatsKMax,
                     PipiIO::LambdaFileName,
                     PipiIO::LambdaHistoName);

    TF1* fit = FitFixedAlphaOnHistogram(hRebinned,
                                        Form("fit_bin8MeV_alpha%.2f", alpha),
                                        alpha,
                                        sourceSizeInit,
                                        PipiIO::FitMin,
                                        PipiIO::FitMax,
                                        PipiIO::BaselineFitMin,
                                        PipiIO::BaselineFitMax);
    if(!fit) continue;

    const double chi2 = fit->GetChisquare();
    const int ndf = fit->GetNDF();
    const double chi2ndf = ndf>0 ? chi2/double(ndf) : 0.0;
    summary.Fill(alpha,
                 fit->GetParameter(0),
                 chi2,
                 ndf,
                 chi2ndf,
                 fit->GetParameter(2),
                 fit->GetParameter(3));

    outputFile->cd();
    fit->Write("", TObject::kOverwrite);
    summary.Write("", TObject::kOverwrite);
    outputFile->Save();
    delete fit;
  }

  outputFile->cd();
  summary.Write("", TObject::kOverwrite);
  outputFile->Close();
  inputFile->Close();
  delete outputFile;
  delete inputFile;
  fracWeakDecay = savedWeakDecay;

  std::cout << "Saved 8 MeV alpha scan to " << outputFileName << "\n";
}

void ScanWeakDecay8MeV(bool Identical,
                       bool WithCoulomb,
                       double fixedAlpha,
                       double weakMin,
                       double weakMax,
                       double weakStep,
                       double fitMin,
                       double fitMax,
                       double baselineFitMin,
                       double baselineFitMax,
                       unsigned catsMomBins)
{
  (void)fitMin;
  (void)fitMax;
  (void)baselineFitMin;
  (void)baselineFitMax;
  (void)catsMomBins;
  (void)fixedAlpha;
  fixedAlpha = sourceAlpha;

  const double savedWeakDecay = fracWeakDecay;

  const TString outputFileName = PipiIO::OutputScanWeakDecay8MeV;

  TFile* inputFile = TFile::Open(PipiIO::DataFileName, "READ");
  if(!inputFile || inputFile->IsZombie()){
    std::cout << "ERROR: cannot open input ROOT file: " << PipiIO::DataFileName << "\n";
    fracWeakDecay = savedWeakDecay;
    return;
  }

  TH1* hDataRaw = dynamic_cast<TH1*>(inputFile->Get(PipiIO::DataHistoName));
  if(!hDataRaw){
    std::cout << "ERROR: cannot find histogram: " << PipiIO::DataHistoName << "\n";
    inputFile->Close();
    fracWeakDecay = savedWeakDecay;
    return;
  }

  TH1D* hRebinned = RebinCfWeighted(hDataRaw, "data_bin8MeV", UniformEdges(0.0, 200.0, 8.0));
  if(!hRebinned){
    inputFile->Close();
    fracWeakDecay = savedWeakDecay;
    return;
  }

  TFile* outputFile = new TFile(outputFileName, "RECREATE");
  TNtuple summary("weakDecay8MeV_scan",
                  "weakDecay8MeV_scan",
                  "weakDecay:alpha:R:chi2:ndf:chi2ndf:p0:p1:lambdaScale");

  outputFile->cd();
  hRebinned->Write();

  const double sourceSizeInit = 12.3;
  for(double weak=weakMin; weak<=weakMax+1e-9; weak+=weakStep){
    fracWeakDecay = weak;
    const double lambdaScale = (1.0-weak)*(1.0-weak);
    std::cout << "ScanWeakDecay8MeV: weak=" << weak
              << ", alpha=" << fixedAlpha << "\n";

    Init_pipi_Theory("CoreReso",
                     Identical,
                     WithCoulomb,
                     0,
                     sourceSizeInit,
                     fixedAlpha,
                     PipiIO::CatsNumMomBins,
                     PipiIO::CatsKMin,
                     PipiIO::CatsKMax,
                     PipiIO::LambdaFileName,
                     PipiIO::LambdaHistoName);

    TF1* fit = FitFixedAlphaOnHistogram(hRebinned,
                                        Form("fit_bin8MeV_weak%.3f_alpha%.2f", weak, fixedAlpha),
                                        fixedAlpha,
                                        sourceSizeInit,
                                        PipiIO::FitMin,
                                        PipiIO::FitMax,
                                        PipiIO::BaselineFitMin,
                                        PipiIO::BaselineFitMax);
    if(!fit) continue;

    const double chi2 = fit->GetChisquare();
    const int ndf = fit->GetNDF();
    const double chi2ndf = ndf>0 ? chi2/double(ndf) : 0.0;
    summary.Fill(weak,
                 fixedAlpha,
                 fit->GetParameter(0),
                 chi2,
                 ndf,
                 chi2ndf,
                 fit->GetParameter(2),
                 fit->GetParameter(3),
                 lambdaScale);

    outputFile->cd();
    fit->Write("", TObject::kOverwrite);
    summary.Write("", TObject::kOverwrite);
    outputFile->Save();
    delete fit;
  }

  outputFile->cd();
  summary.Write("", TObject::kOverwrite);
  outputFile->Close();
  inputFile->Close();
  delete outputFile;
  delete inputFile;
  fracWeakDecay = savedWeakDecay;

  std::cout << "Saved 8 MeV weak-decay scan to " << outputFileName << "\n";
}

// chi2_info: [0]=chi2, [1]=ndf, [2]=chi2/ndf
// AfterFitPara: [0..4]=baseline params, [5]=femto param (or lambda/fhalo)
// baselineFunc: double (*)(double*, double*)
// femtoFunc:    double (*)(double*, double*)
// =========================================================
// Fit driver
// =========================================================
void fitData3DTo1D(bool Identical,
                   bool WithCoulomb,
                   bool usePRF3DShift)
{

  const TString outputFileName = PipiIO::OutputFit3DTo1D;

  TFile* FileROOT = TFile::Open(PipiIO::DataFileName,"READ");
  if(!FileROOT || FileROOT->IsZombie()){
    std::cout << "ERROR: cannot open input ROOT file: " << PipiIO::DataFileName << "\n";
    return;
  }

  TH1* hData = dynamic_cast<TH1*>(FileROOT->Get(PipiIO::DataHistoName));
  if(!hData){
    std::cout << "ERROR: cannot find histogram: " << PipiIO::DataHistoName << "\n";
    FileROOT->Close();
    return;
  }

  TH1* hDataToFit = (TH1*)hData->Clone("copyhisto_3DTo1D");
  if(!hDataToFit){
    std::cout << "ERROR: cannot clone input histogram\n";
    FileROOT->Close();
    return;
  }

  const double RoutInit = 6.0;
  const double RsideInit = 6.0;
  const double RlongInit = 6.0;
  const double SourceAlphaInit = sourceAlpha;

  const TString sourceType = usePRF3DShift ? "CoreReso3DTo1D_PRF" : "CoreReso3DTo1D";

  Init_pipi_Theory(sourceType,
                   Identical,
                   WithCoulomb,
                   0,
                   RoutInit,
                   SourceAlphaInit,
                   PipiIO::CatsNumMomBins,
                   PipiIO::CatsKMin,
                   PipiIO::CatsKMax,
                   PipiIO::LambdaFileName, PipiIO::LambdaHistoName);
  if(!gFemto_pipi_fitter){
    std::cout << "ERROR: Init_pipi_Theory failed, skip 3D-to-1D fit.\n";
    FileROOT->Close();
    return;
  }

  const double fitMin = PipiIO::FitMin;
  const double fitMax = PipiIO::FitMax;
  TF1* fit_pipi = new TF1("fit_pipi_3DTo1D", MickeyFitter_pipi_3DTo1D_Pol1, fitMin, fitMax, 6);
  fit_pipi->SetNpx(PipiIO::FitNpx);

  const double baselineFitMin = PipiIO::BaselineFitMin;
  const double baselineFitMax = PipiIO::BaselineFitMax;
  TF1* baseline_prefit = new TF1("baseline_prefit_3DTo1D", "[0]*(1.0+[1]*x)", baselineFitMin, baselineFitMax);
  baseline_prefit->SetParameter(0, 1.0);
  baseline_prefit->SetParameter(1, 0.0);
  hDataToFit->Fit(baseline_prefit, "RQ0", "", baselineFitMin, baselineFitMax);

  fit_pipi->SetParameter(0, RoutInit);
  fit_pipi->SetParLimits(0, 5.0, 15.0);
  fit_pipi->SetParameter(1, RsideInit);
  fit_pipi->SetParLimits(1, 5.0, 15.0);
  fit_pipi->SetParameter(2, RlongInit);
  fit_pipi->SetParLimits(2, 5.0, 15.0);

  fit_pipi->SetParameter(3, SourceAlphaInit);
  fit_pipi->FixParameter(3, SourceAlphaInit);

  fit_pipi->SetParameter(4, baseline_prefit->GetParameter(0));
  fit_pipi->FixParameter(4, baseline_prefit->GetParameter(0));
  fit_pipi->SetParameter(5, baseline_prefit->GetParameter(1));
  fit_pipi->FixParameter(5, baseline_prefit->GetParameter(1));

  TGraphAsymmErrors* grBand = new TGraphAsymmErrors();
  TGraphErrors* grfit = new TGraphErrors();
  TGraph gUp;
  TGraph gDn;
  TGraph gCtrl;
  std::cout << "start 3D-to-1D Fit !!!\n";
  DrawErrorBand(hDataToFit, fit_pipi, grfit, grBand, gUp, gDn, gCtrl, 0.95);

  std::cout << "3D-to-1D Fit done:\n";
  std::cout << "SourceType = " << sourceType << "\n";
  std::cout << "Rout = " << fit_pipi->GetParameter(0) << "\n";
  std::cout << "Rside = " << fit_pipi->GetParameter(1) << "\n";
  std::cout << "Rlong = " << fit_pipi->GetParameter(2) << "\n";
  std::cout << "alpha = " << fit_pipi->GetParameter(3) << "\n";
  std::cout << "fracWeakDecay = " << fracWeakDecay << "\n";
  std::cout << "Fit range = [" << fitMin << ", " << fitMax << "] MeV/c\n";
  std::cout << "Baseline p0 fixed from [" << baselineFitMin << ", " << baselineFitMax
            << "] = " << fit_pipi->GetParameter(4) << "\n";
  std::cout << "Baseline p1 fixed from [" << baselineFitMin << ", " << baselineFitMax
            << "] = " << fit_pipi->GetParameter(5) << "\n";

  TFile* SaveFile = new TFile(outputFileName, "RECREATE");
  SaveFile->cd();
  hDataToFit->Write();
  baseline_prefit->Write();
  fit_pipi->SetName("CATS_fit_3DTo1D");
  fit_pipi->Write();

  TF1* fit_pipi_baseline = new TF1("fit_pipi_baseline_3DTo1D", baseline3, 0, 400, 4);
  fit_pipi_baseline->SetParameter(0, fit_pipi->GetParameter(4));
  fit_pipi_baseline->SetParameter(1, fit_pipi->GetParameter(5));
  fit_pipi_baseline->SetParameter(2, 0.0);
  fit_pipi_baseline->SetParameter(3, 0.0);
  fit_pipi_baseline->SetLineColor(kGreen-2);
  fit_pipi_baseline->Write();

  SaveFile->Close();
  FileROOT->Close();

  std::cout << "Saved 3D-to-1D fit to " << outputFileName << "\n";

  delete grBand;
  delete grfit;
  delete baseline_prefit;
  delete fit_pipi_baseline;
  delete SaveFile;
  delete FileROOT;
}

void fitData(bool Identical, bool WithCoulomb)
{

  const TString outputFileName = PipiIO::OutputFit1D;

  TFile* FileROOT = TFile::Open(PipiIO::DataFileName,"READ");
  if(!FileROOT || FileROOT->IsZombie()){
    std::cout << "ERROR: cannot open input ROOT file: " << PipiIO::DataFileName << "\n";
    return;
  }

  TH1* hData = dynamic_cast<TH1*>(FileROOT->Get(PipiIO::DataHistoName));
  if(!hData){
    std::cout << "ERROR: cannot find histogram: " << PipiIO::DataHistoName << "\n";
    return;
  }

  TH1* hDataToFit = (TH1*)hData->Clone("copyhisto");
  if(!hDataToFit){
    std::cout << "ERROR: cannot clone input histogram\n";
    return;
  }

  const double SourceSizeInit = 12;
  const double SourceAlphaInit = sourceAlpha;
  Init_pipi_Theory("CoreResoDoubleGauss",
                   Identical,
                   WithCoulomb,
                   0,
                   SourceSizeInit,
                   SourceAlphaInit,
                   PipiIO::CatsNumMomBins,
                   PipiIO::CatsKMin,
                   PipiIO::CatsKMax,
                   PipiIO::LambdaFileName, PipiIO::LambdaHistoName);
  if(!gFemto_pipi_fitter){
    std::cout << "ERROR: Init_pipi_Theory failed, skip 1D fit.\n";
    FileROOT->Close();
    return;
  }

  const double fitMin = PipiIO::FitMin;
  const double fitMax = PipiIO::FitMax;
  TF1* fit_pipi = new TF1("fit_pipi", MickeyFitter_pipi_Pol3, fitMin, fitMax, 4);
  fit_pipi->SetNpx(PipiIO::FitNpx);

  // Determine the non-femtoscopic baseline away from the low-k femtoscopic
  // signal, then keep it fixed in the CATS fit.
  const double baselineFitMin = PipiIO::BaselineFitMin;
  const double baselineFitMax = PipiIO::BaselineFitMax;
  TF1* baseline_prefit = new TF1("baseline_prefit", "[0]*(1.0+[1]*x)", baselineFitMin, baselineFitMax);
  baseline_prefit->SetParameter(0, 1.0);
  baseline_prefit->SetParameter(1, 0.0);
  hDataToFit->Fit(baseline_prefit, "RQ0", "", baselineFitMin, baselineFitMax);

  fit_pipi->SetParameter(0, SourceSizeInit);
  fit_pipi->SetParLimits(0, 3.0, 18.0);

  fit_pipi->SetParameter(1, SourceAlphaInit);
  fit_pipi->FixParameter(1, SourceAlphaInit);

  fit_pipi->SetParameter(2, baseline_prefit->GetParameter(0));
  fit_pipi->FixParameter(2, baseline_prefit->GetParameter(0));
  fit_pipi->SetParameter(3, baseline_prefit->GetParameter(1));
  fit_pipi->FixParameter(3, baseline_prefit->GetParameter(1));

  TGraphAsymmErrors* grBand = new TGraphAsymmErrors();
  TGraphErrors* grfit = new TGraphErrors();
  TGraph gUp;
  TGraph gDn;
  TGraph gCtrl;
  std::cout << "start Fit !!!\n";
  DrawErrorBand(hDataToFit, fit_pipi, grfit, grBand, gUp, gDn, gCtrl, 0.95);

  std::cout << "Fit done:\n";
  std::cout << "R = " << fit_pipi->GetParameter(0) << "\n";
  std::cout << "double Gaussian tail width = " << doubleGaussTailWidth << " fm (fixed)\n";
  std::cout << "double Gaussian core fraction = " << doubleGaussCoreFraction << " (fixed)\n";
  std::cout << "fracWeakDecay = " << fracWeakDecay << "\n";
  std::cout << "Fit range = [" << fitMin << ", " << fitMax << "] MeV/c\n";
  std::cout << "Baseline p0 fixed from [" << baselineFitMin << ", " << baselineFitMax
            << "] = " << fit_pipi->GetParameter(2) << "\n";
  std::cout << "Baseline p1 fixed from [" << baselineFitMin << ", " << baselineFitMax
            << "] = " << fit_pipi->GetParameter(3) << "\n";

  TFile* SaveFile = new TFile(outputFileName, "RECREATE");
  SaveFile->cd();

  hDataToFit->Write();

  fit_pipi->SetName("CATS_fit");
  fit_pipi->Write();

  TF1* fit_pipi_baseline = new TF1("fit_pipi_baseline", baseline3, 0, 200, 4);
  fit_pipi_baseline->SetParameter(0, fit_pipi->GetParameter(2));
  fit_pipi_baseline->SetParameter(1, fit_pipi->GetParameter(3));
  fit_pipi_baseline->SetParameter(2, 0.0);
  fit_pipi_baseline->SetParameter(3, 0.0);
  SaveFile->cd();
  fit_pipi_baseline->SetLineColor(kGreen-2);
  fit_pipi_baseline->Write();

  baseline_prefit->Write();

  printf("Saved: %s\n", outputFileName.Data());

  SaveFile->Close();

}
