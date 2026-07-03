#include "CATSCk.h"

#include "Basics.h"

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
#include "TH1D.h"
#include "TH2F.h"
#include "TLegend.h"
#include "TNtuple.h"
#include "TSystem.h"

#include <cmath>
#include <iostream>
#include <vector>

namespace {

enum class FlowMode {
  None,
  V2YieldOnly,
  V2YieldAndGeometry
};

void SetGaussianSource(CATS& kitty, const double sourceSize)
{
  CATSparameters sourcePars(CATSparameters::tSource, 1, true);
  sourcePars.SetParameter(0, sourceSize);
  kitty.SetAnaSource(GaussSource, sourcePars);
  kitty.SetUseAnalyticSource(true);
  kitty.SetAutoNormSource(false);
}

double PositiveFlowScale(const double phi, const double sourceAnisotropy)
{
  const double scale = 1.0 + 2.0 * sourceAnisotropy * std::cos(2.0 * phi);
  return scale > 0.05 ? scale : 0.05;
}

bool AcceptByV2(DLM_Random& ranGen, const double phi, const double v2)
{
  const double weight = 1.0 + 2.0 * v2 * std::cos(2.0 * phi);
  const double maxWeight = 1.0 + 2.0 * std::abs(v2);
  return ranGen.Uniform(0.0, maxWeight) <= weight;
}

DLM_CleverMcLevyResoTM* SetResonanceSource(CATS& kitty,
                                           const double sourceSize,
                                           const double alpha,
                                           const double shortLivedResonanceFraction,
                                           const double resonanceTau,
                                           const FlowMode flowMode,
                                           const double v2,
                                           const double sourceAnisotropy,
                                           const double cutoff,
                                           const TString& primaryResonanceFile,
                                           const TString& resonanceResonanceFile)
{
  DLM_CleverMcLevyResoTM* source = new DLM_CleverMcLevyResoTM();
  source->InitStability(1, alpha - 1e-6, alpha + 1e-6);
  source->InitScale(160, 4.0, 30.0);
  source->InitRad(400, 1e-4, 128);
  source->InitType(2);

  source->SetUpReso(0, shortLivedResonanceFraction);
  source->SetUpReso(1, shortLivedResonanceFraction);

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

  TFile* filePR = TFile::Open(primaryResonanceFile, "READ");
  if(!filePR || filePR->IsZombie()) {
    std::cout << "ERROR: cannot open primary-resonance file: "
              << primaryResonanceFile << "\n";
    if(filePR) {
      filePR->Close();
      delete filePR;
    }
  }
  else {
    TNtuple* tuplePR = dynamic_cast<TNtuple*>(filePR->Get("InfoTuple_ClosePairs"));
    if(!tuplePR) {
      std::cout << "ERROR: cannot find InfoTuple_ClosePairs in "
                << primaryResonanceFile << "\n";
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

        const double phi = ranGen.Uniform(0.0, 2.0 * M_PI);
        if(flowMode != FlowMode::None && !AcceptByV2(ranGen, phi, v2)) continue;

        const double flowScale = flowMode == FlowMode::V2YieldAndGeometry
                               ? PositiveFlowScale(phi, sourceAnisotropy)
                               : 1.0;
        // const double shift = flowScale * ranGen.Exponential(fM2 / (fP2 * resonanceTau));
        const double shift = ranGen.Exponential(fM2 / (fP2 * resonanceTau));
        source->AddBGT_PR(shift, -std::cos(angleRcP2));
        source->AddBGT_RP(shift,  std::cos(angleRcP2));
        ++addedPR;
      }
      std::cout << "Loaded " << addedPR << " PR/RP resonance kernels from "
                << primaryResonanceFile << "\n";
    }
    filePR->Close();
    delete filePR;
  }

  TFile* fileRR = TFile::Open(resonanceResonanceFile, "READ");
  if(!fileRR || fileRR->IsZombie()) {
    std::cout << "ERROR: cannot open resonance-resonance file: "
              << resonanceResonanceFile << "\n";
    if(fileRR) {
      fileRR->Close();
      delete fileRR;
    }
  }
  else {
    TNtuple* tupleRR = dynamic_cast<TNtuple*>(fileRR->Get("InfoTuple_ClosePairs"));
    if(!tupleRR) {
      std::cout << "ERROR: cannot find InfoTuple_ClosePairs in "
                << resonanceResonanceFile << "\n";
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

        const double phi1 = ranGen.Uniform(0.0, 2.0 * M_PI);
        const double phi2 = ranGen.Uniform(0.0, 2.0 * M_PI);
        if(flowMode != FlowMode::None &&
           (!AcceptByV2(ranGen, phi1, v2) || !AcceptByV2(ranGen, phi2, v2))) {
          continue;
        }

        const double flowScale1 = flowMode == FlowMode::V2YieldAndGeometry
                                ? PositiveFlowScale(phi1, sourceAnisotropy)
                                : 1.0;
        const double flowScale2 = flowMode == FlowMode::V2YieldAndGeometry
                                ? PositiveFlowScale(phi2, sourceAnisotropy)
                                : 1.0;
        // const double shift1 = flowScale1 * ranGen.Exponential(fM1 / (fP1 * resonanceTau));
        // const double shift2 = flowScale2 * ranGen.Exponential(fM2 / (fP2 * resonanceTau));

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
                << resonanceResonanceFile << "\n";
    }
    fileRR->Close();
    delete fileRR;
  }

  if(alpha == 1.0 || alpha == 2.0) {
    source->InitNumMcIter(1000000);
  }
  else {
    source->InitNumMcIter(100000);
  }

  kitty.SetAnaSource(CatsSourceForwarder, source, 2);
  kitty.SetAnaSource(0, sourceSize);
  kitty.SetAnaSource(1, alpha);
  kitty.SetUseAnalyticSource(true);
  kitty.SetAutoNormSource(false);

  return source;
}

template <class SourceType>
void LoadOneDimensionalResonanceKernels(SourceType* source,
                                        const double resonanceTau,
                                        const double cutoff,
                                        const FlowMode flowMode,
                                        const double v2,
                                        const double sourceAnisotropy,
                                        const TString& primaryResonanceFile,
                                        const TString& resonanceResonanceFile)
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

  TFile* filePR = TFile::Open(primaryResonanceFile, "READ");
  if(!filePR || filePR->IsZombie()) {
    std::cout << "ERROR: cannot open primary-resonance file: "
              << primaryResonanceFile << "\n";
    if(filePR) {
      filePR->Close();
      delete filePR;
    }
  }
  else {
    TNtuple* tuplePR = dynamic_cast<TNtuple*>(filePR->Get("InfoTuple_ClosePairs"));
    if(!tuplePR) {
      std::cout << "ERROR: cannot find InfoTuple_ClosePairs in "
                << primaryResonanceFile << "\n";
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

        const double phi = ranGen.Uniform(0.0, 2.0 * M_PI);
        if(flowMode != FlowMode::None && !AcceptByV2(ranGen, phi, v2)) continue;
        const double flowScale = flowMode == FlowMode::V2YieldAndGeometry
                               ? PositiveFlowScale(phi, sourceAnisotropy)
                               : 1.0;
        (void)flowScale;

        const double shift = ranGen.Exponential(fM2 / (fP2 * resonanceTau));
        source->AddBGT_PR(shift, -std::cos(angleRcP2));
        source->AddBGT_RP(shift,  std::cos(angleRcP2));
        ++addedPR;
      }
      std::cout << "Loaded " << addedPR << " PR/RP resonance kernels from "
                << primaryResonanceFile << "\n";
    }
    filePR->Close();
    delete filePR;
  }

  TFile* fileRR = TFile::Open(resonanceResonanceFile, "READ");
  if(!fileRR || fileRR->IsZombie()) {
    std::cout << "ERROR: cannot open resonance-resonance file: "
              << resonanceResonanceFile << "\n";
    if(fileRR) {
      fileRR->Close();
      delete fileRR;
    }
  }
  else {
    TNtuple* tupleRR = dynamic_cast<TNtuple*>(fileRR->Get("InfoTuple_ClosePairs"));
    if(!tupleRR) {
      std::cout << "ERROR: cannot find InfoTuple_ClosePairs in "
                << resonanceResonanceFile << "\n";
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

        const double phi1 = ranGen.Uniform(0.0, 2.0 * M_PI);
        const double phi2 = ranGen.Uniform(0.0, 2.0 * M_PI);
        if(flowMode != FlowMode::None &&
           (!AcceptByV2(ranGen, phi1, v2) || !AcceptByV2(ranGen, phi2, v2))) {
          continue;
        }
        const double flowScale1 = flowMode == FlowMode::V2YieldAndGeometry
                                ? PositiveFlowScale(phi1, sourceAnisotropy)
                                : 1.0;
        const double flowScale2 = flowMode == FlowMode::V2YieldAndGeometry
                                ? PositiveFlowScale(phi2, sourceAnisotropy)
                                : 1.0;
        (void)flowScale1;
        (void)flowScale2;

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
                << resonanceResonanceFile << "\n";
    }
    fileRR->Close();
    delete fileRR;
  }
}

DLM_CleverMcNonGaussResoTM* SetDoubleGaussianResonanceSource(
    CATS& kitty,
    const double sourceCoreSize,
    const double sourceTailSize,
    const double sourceCoreFraction,
    const double shortLivedResonanceFraction,
    const double resonanceTau,
    const double cutoff,
    const TString& primaryResonanceFile,
    const TString& resonanceResonanceFile)
{
  DLM_CleverMcNonGaussResoTM* source = new DLM_CleverMcNonGaussResoTM();
  source->InitScale(120, 2.0, 30.0);
  source->InitTailScale(120, 5.0, 160.0);
  source->InitCoreFraction(80, 0.0, 1.0);
  source->InitRad(500, 1e-4, 220.0);
  source->InitType(1);

  source->SetUpReso(0, shortLivedResonanceFraction);
  source->SetUpReso(1, shortLivedResonanceFraction);

  LoadOneDimensionalResonanceKernels(source,
                                     resonanceTau,
                                     cutoff,
                                     FlowMode::None,
                                     0.0,
                                     0.0,
                                     primaryResonanceFile,
                                     resonanceResonanceFile);
  source->InitNumMcIter(1000000);

  kitty.SetAnaSource(CatsSourceForwarder, source, 3);
  kitty.SetAnaSource(0, sourceCoreSize);
  kitty.SetAnaSource(1, sourceTailSize);
  kitty.SetAnaSource(2, sourceCoreFraction);
  kitty.SetUseAnalyticSource(true);
  kitty.SetAutoNormSource(false);

  return source;
}

void SetPipiSystem(CATS& kitty,
                   const unsigned numMomBins,
                   const double kMin,
                   const double kMax,
                   const bool identical,
                   const bool withCoulomb)
{
  kitty.SetMomBins(numMomBins, kMin, kMax);

  kitty.SetNumChannels(1);
  kitty.SetNumPW(0, 0);
  kitty.SetSpin(0, 0);
  kitty.SetChannelWeight(0, 1.0);

  kitty.SetQ1Q2(withCoulomb ? 1 : 0);
  kitty.SetPdgId(identical ? 211 : -211, 211);

  if(withCoulomb) {
    kitty.SetRedMass(0.5 * Mass_Pic);
  }
  else if(identical) {
    kitty.SetRedMass(0.5 * Mass_Pi0);
  }
  else {
    kitty.SetRedMass((Mass_Pi0 * Mass_Pic) / (Mass_Pi0 + Mass_Pic));
  }
}

TH1D* MakeHistogramFromCk(const TString& histName, DLM_Ck& ck)
{
  const unsigned numBins = ck.GetNbins();
  if(numBins < 1) {
    std::cout << "ERROR: CATS returned an empty C(k*) object.\n";
    return nullptr;
  }

  const double firstCenter = ck.GetBinCenter(0, 0);
  double binWidth = 1.0;
  if(numBins > 1) {
    binWidth = ck.GetBinCenter(0, 1) - firstCenter;
  }

  const double xMin = firstCenter - 0.5 * binWidth;
  const double xMax = xMin + numBins * binWidth;

  TH1D* hist = new TH1D(histName, histName, numBins, xMin, xMax);
  hist->SetDirectory(nullptr);
  hist->SetStats(false);
  hist->GetXaxis()->SetTitle("k* (MeV/c)");
  hist->GetYaxis()->SetTitle("C(k*)");
  hist->SetLineColor(kRed + 1);
  hist->SetLineWidth(2);

  for(unsigned iBin = 0; iBin < numBins; ++iBin) {
    hist->SetBinContent(iBin + 1, ck.GetBinContent(iBin));
  }

  return hist;
}

DLM_Histo<double>* MakeFeedDownFractionFromLambda(const TH1* hLambda,
                                                  const double weakDecayFraction,
                                                  const double xScale)
{
  if(!hLambda) return nullptr;

  const unsigned numBins = hLambda->GetNbinsX();
  std::vector<double> edges(numBins + 1);
  for(unsigned iBin = 0; iBin <= numBins; ++iBin) {
    edges[iBin] = hLambda->GetXaxis()->GetBinLowEdge(iBin + 1) * xScale;
  }

  DLM_Histo<double>* hFeedDown = new DLM_Histo<double>();
  hFeedDown->SetUp(1);
  hFeedDown->SetUp(0, numBins, edges.data());
  hFeedDown->Initialize(true);

  unsigned whichBin[1];
  const double weakSurvival = (1.0 - weakDecayFraction) * (1.0 - weakDecayFraction);

  for(unsigned iBin = 0; iBin < numBins; ++iBin) {
    const double lambdaMain = hLambda->GetBinContent(iBin + 1) * weakSurvival;
    const double lambdaErr = hLambda->GetBinError(iBin + 1) * weakSurvival;

    if(!std::isfinite(lambdaMain) || lambdaMain < 0.0 || lambdaMain > 1.0) {
      std::cout << "ERROR: invalid lambda value in bin "
                << iBin + 1 << ": " << lambdaMain << "\n";
      delete hFeedDown;
      return nullptr;
    }

    whichBin[0] = iBin;
    hFeedDown->SetBinContent(whichBin, 1.0 - lambdaMain);
    hFeedDown->SetBinError(whichBin, lambdaErr);
  }

  return hFeedDown;
}

DLM_Histo<double>* LoadFeedDownFraction(const TString& fileName,
                                        const TString& histName,
                                        const double weakDecayFraction,
                                        const double xScale)
{
  TFile* file = TFile::Open(fileName, "READ");
  if(!file || file->IsZombie()) {
    std::cout << "ERROR: cannot open lambda ROOT file: "
              << fileName << "\n";
    if(file) {
      file->Close();
      delete file;
    }
    return nullptr;
  }

  TH1* hLambda = dynamic_cast<TH1*>(file->Get(histName));
  if(!hLambda) {
    std::cout << "ERROR: cannot find lambda histogram "
              << histName << " in " << fileName << "\n";
    file->Close();
    delete file;
    return nullptr;
  }

  DLM_Histo<double>* hFeedDown =
      MakeFeedDownFractionFromLambda(hLambda, weakDecayFraction, xScale);

  file->Close();
  delete file;
  return hFeedDown;
}

TH1D* MakeHistogramFromDecomposition(const TString& histName,
                                     DLM_CkDecomposition& ckDecomposition,
                                     DLM_Ck& ckTemplate)
{
  TH1D* hist = MakeHistogramFromCk(histName, ckTemplate);
  if(!hist) return nullptr;

  for(int iBin = 1; iBin <= hist->GetNbinsX(); ++iBin) {
    hist->SetBinContent(iBin, ckDecomposition.EvalCk(hist->GetBinCenter(iBin)));
  }

  hist->SetLineColor(kRed + 1);
  hist->SetLineWidth(2);
  return hist;
}

TH1D* BuildResonanceCurve(const TString& histName,
                          const bool identical,
                          const bool withCoulomb,
                          const double sourceSize,
                          const unsigned numMomBins,
                          const double kMin,
                          const double kMax,
                          const TString& lambdaFileName,
                          const TString& lambdaHistName,
                          const double weakDecayFraction,
                          const double lambdaXScale,
                          const double sourceAlpha,
                          const double shortLivedResonanceFraction,
                          const double resonanceTau,
                          const FlowMode flowMode,
                          const double v2,
                          const double sourceAnisotropy,
                          const TString& primaryResonanceFile,
                          const TString& resonanceResonanceFile)
{
  CATS kitty;
  DLM_CleverMcLevyResoTM* resonanceSource =
      SetResonanceSource(kitty,
                         sourceSize,
                         sourceAlpha,
                         shortLivedResonanceFraction,
                         resonanceTau,
                         flowMode,
                         v2,
                         sourceAnisotropy,
                         kMax,
                         primaryResonanceFile,
                         resonanceResonanceFile);

  SetPipiSystem(kitty, numMomBins, kMin, kMax, identical, withCoulomb);
  kitty.SetNotifications(CATS::nWarning);
  kitty.SetNotifications(CATS::nError);
  kitty.KillTheCat();

  DLM_Ck ck(kitty.GetNumSourcePars(), 0, kitty);
  ck.Update(true);

  TH1D* hist = nullptr;
  if(lambdaFileName.Length() > 0 && lambdaHistName.Length() > 0) {
    DLM_Histo<double>* hFeedDown =
        LoadFeedDownFraction(lambdaFileName,
                             lambdaHistName,
                             weakDecayFraction,
                             lambdaXScale);
    if(hFeedDown) {
      DLM_CkDecomposition ckDecomposition(histName, 1, ck, nullptr);
      ckDecomposition.AddContribution(0, *hFeedDown, DLM_CkDecomposition::cFeedDown);
      ckDecomposition.Update(true, true);
      hist = MakeHistogramFromDecomposition(histName, ckDecomposition, ck);
      delete hFeedDown;
    }
  }
  else {
    hist = MakeHistogramFromCk(histName, ck);
  }

  delete resonanceSource;
  return hist;
}

} // namespace

void DrawCATSCk(bool identical,
                bool withCoulomb,
                double sourceSize,
                unsigned numMomBins,
                double kMin,
                double kMax,
                TString outputFileName,
                TString histName,
                TString lambdaFileName,
                TString lambdaHistName,
                double weakDecayFraction,
                double lambdaXScale,
                bool useResonanceSource,
                double sourceAlpha,
                double shortLivedResonanceFraction,
                double resonanceTau,
                TString primaryResonanceFile,
                TString resonanceResonanceFile)
{
  gSystem->mkdir(gSystem->DirName(outputFileName), true);

  CATS kitty;
  DLM_CleverMcLevyResoTM* resonanceSource = nullptr;
  if(useResonanceSource) {
    resonanceSource = SetResonanceSource(kitty,
                                         sourceSize,
                                         sourceAlpha,
                                         shortLivedResonanceFraction,
                                         resonanceTau,
                                         FlowMode::None,
                                         0.0,
                                         0.0,
                                         kMax,
                                         primaryResonanceFile,
                                         resonanceResonanceFile);
  }
  else {
    SetGaussianSource(kitty, sourceSize);
  }

  SetPipiSystem(kitty, numMomBins, kMin, kMax, identical, withCoulomb);
  kitty.SetNotifications(CATS::nWarning);
  kitty.SetNotifications(CATS::nError);
  kitty.KillTheCat();

  DLM_Ck ck(kitty.GetNumSourcePars(), 0, kitty);
  ck.Update(true);

  TH1D* hCore = MakeHistogramFromCk(histName + "_core", ck);
  if(!hCore) return;

  TH1D* hist = nullptr;
  DLM_Histo<double>* hFeedDown = nullptr;

  if(lambdaFileName.Length() > 0 && lambdaHistName.Length() > 0) {
    hFeedDown = LoadFeedDownFraction(lambdaFileName,
                                     lambdaHistName,
                                     weakDecayFraction,
                                     lambdaXScale);
    if(!hFeedDown) {
      delete hCore;
      return;
    }

    DLM_CkDecomposition ckDecomposition("pipi_lambda_feeddown", 1, ck, nullptr);
    ckDecomposition.AddContribution(0, *hFeedDown, DLM_CkDecomposition::cFeedDown);
    ckDecomposition.Update(true, true);
    hist = MakeHistogramFromDecomposition(histName, ckDecomposition, ck);
  }
  else {
    hist = static_cast<TH1D*>(hCore->Clone(histName));
    hist->SetDirectory(nullptr);
  }

  if(!hist) {
    delete hFeedDown;
    delete hCore;
    return;
  }

  TCanvas canvas("cCATS_Ck", "CATS C(k*)", 900, 650);
  hCore->SetLineColor(kGray + 2);
  hCore->SetLineStyle(2);
  hCore->SetLineWidth(2);

  if(hFeedDown) {
    hCore->Draw("hist l");
    hist->Draw("hist l same");

    TLegend legend(0.58, 0.72, 0.88, 0.88);
    legend.SetBorderSize(0);
    legend.SetFillStyle(0);
    legend.AddEntry(hCore, "CATS core", "l");
    legend.AddEntry(hist, "CATS + lambda feed-down", "l");
    legend.Draw();
  }
  else {
    hist->Draw("hist l");
  }

  TFile outputFile(outputFileName, "RECREATE");
  hCore->Write("", TObject::kOverwrite);
  hist->Write("", TObject::kOverwrite);
  canvas.Write("", TObject::kOverwrite);
  outputFile.Close();

  TString pdfName = outputFileName;
  TString pngName = outputFileName;
  canvas.SaveAs(pdfName.ReplaceAll(".root", ".pdf"));
  canvas.SaveAs(pngName.ReplaceAll(".root", ".png"));

  std::cout << "Saved CATS correlation function";
  if(resonanceSource) {
    std::cout << " with resonance source reshaping";
  }
  if(hFeedDown) {
    std::cout << " with lambda feed-down";
  }
  std::cout << " to " << outputFileName << "\n";

  delete hFeedDown;
  delete hCore;
  delete hist;
  delete resonanceSource;
}

void DrawCATSCkDoubleGaussianReso(bool identical,
                                  bool withCoulomb,
                                  double sourceCoreSize,
                                  double sourceTailSize,
                                  double sourceCoreFraction,
                                  unsigned numMomBins,
                                  double kMin,
                                  double kMax,
                                  TString outputFileName,
                                  TString histName,
                                  TString lambdaFileName,
                                  TString lambdaHistName,
                                  double weakDecayFraction,
                                  double lambdaXScale,
                                  double shortLivedResonanceFraction,
                                  double resonanceTau,
                                  TString primaryResonanceFile,
                                  TString resonanceResonanceFile)
{
  gSystem->mkdir(gSystem->DirName(outputFileName), true);

  CATS kitty;
  DLM_CleverMcNonGaussResoTM* resonanceSource =
      SetDoubleGaussianResonanceSource(kitty,
                                       sourceCoreSize,
                                       sourceTailSize,
                                       sourceCoreFraction,
                                       shortLivedResonanceFraction,
                                       resonanceTau,
                                       kMax,
                                       primaryResonanceFile,
                                       resonanceResonanceFile);

  SetPipiSystem(kitty, numMomBins, kMin, kMax, identical, withCoulomb);
  kitty.SetNotifications(CATS::nWarning);
  kitty.SetNotifications(CATS::nError);
  kitty.KillTheCat();

  DLM_Ck ck(kitty.GetNumSourcePars(), 0, kitty);
  ck.Update(true);

  TH1D* hCore = MakeHistogramFromCk(histName + "_core", ck);
  if(!hCore) {
    delete resonanceSource;
    return;
  }

  TH1D* hist = nullptr;
  DLM_Histo<double>* hFeedDown = nullptr;

  if(lambdaFileName.Length() > 0 && lambdaHistName.Length() > 0) {
    hFeedDown = LoadFeedDownFraction(lambdaFileName,
                                     lambdaHistName,
                                     weakDecayFraction,
                                     lambdaXScale);
    if(!hFeedDown) {
      delete hCore;
      delete resonanceSource;
      return;
    }

    DLM_CkDecomposition ckDecomposition("pipi_double_gauss_lambda_feeddown", 1, ck, nullptr);
    ckDecomposition.AddContribution(0, *hFeedDown, DLM_CkDecomposition::cFeedDown);
    ckDecomposition.Update(true, true);
    hist = MakeHistogramFromDecomposition(histName, ckDecomposition, ck);
  }
  else {
    hist = static_cast<TH1D*>(hCore->Clone(histName));
    hist->SetDirectory(nullptr);
  }

  if(!hist) {
    delete hFeedDown;
    delete hCore;
    delete resonanceSource;
    return;
  }

  TCanvas canvas("cCATS_Ck_DoubleGaussian", "CATS C(k*) double Gaussian source", 900, 650);
  hCore->SetLineColor(kGray + 2);
  hCore->SetLineStyle(2);
  hCore->SetLineWidth(2);

  if(hFeedDown) {
    hCore->Draw("hist l");
    hist->Draw("hist l same");

    TLegend legend(0.52, 0.68, 0.88, 0.88);
    legend.SetBorderSize(0);
    legend.SetFillStyle(0);
    legend.AddEntry(hCore, "CATS double-Gaussian core", "l");
    legend.AddEntry(hist, "CATS + lambda feed-down", "l");
    legend.Draw();
  }
  else {
    hist->Draw("hist l");
  }

  TFile outputFile(outputFileName, "RECREATE");
  hCore->Write("", TObject::kOverwrite);
  hist->Write("", TObject::kOverwrite);
  canvas.Write("", TObject::kOverwrite);
  outputFile.Close();

  TString pdfName = outputFileName;
  TString pngName = outputFileName;
  canvas.SaveAs(pdfName.ReplaceAll(".root", ".pdf"));
  canvas.SaveAs(pngName.ReplaceAll(".root", ".png"));

  std::cout << "Saved CATS double-Gaussian resonance source correlation function to "
            << outputFileName << "\n";
  std::cout << "Double-Gaussian source: R_core=" << sourceCoreSize
            << " fm, R_tail=" << sourceTailSize
            << " fm, f_core=" << sourceCoreFraction << "\n";

  delete hFeedDown;
  delete hCore;
  delete hist;
  delete resonanceSource;
}

void CompareFlowModulations(bool identical,
                            bool withCoulomb,
                            double sourceSize,
                            unsigned numMomBins,
                            double kMin,
                            double kMax,
                            TString outputFileName,
                            TString lambdaFileName,
                            TString lambdaHistName,
                            double weakDecayFraction,
                            double lambdaXScale,
                            double sourceAlpha,
                            double shortLivedResonanceFraction,
                            double resonanceTau,
                            double v2,
                            double sourceAnisotropy,
                            TString primaryResonanceFile,
                            TString resonanceResonanceFile)
{
  gSystem->mkdir(gSystem->DirName(outputFileName), true);

  TH1D* hNoFlow = BuildResonanceCurve("hCATS_noFlow",
                                      identical,
                                      withCoulomb,
                                      sourceSize,
                                      numMomBins,
                                      kMin,
                                      kMax,
                                      lambdaFileName,
                                      lambdaHistName,
                                      weakDecayFraction,
                                      lambdaXScale,
                                      sourceAlpha,
                                      shortLivedResonanceFraction,
                                      resonanceTau,
                                      FlowMode::None,
                                      v2,
                                      sourceAnisotropy,
                                      primaryResonanceFile,
                                      resonanceResonanceFile);

  TH1D* hV2Only = BuildResonanceCurve("hCATS_v2YieldOnly",
                                      identical,
                                      withCoulomb,
                                      sourceSize,
                                      numMomBins,
                                      kMin,
                                      kMax,
                                      lambdaFileName,
                                      lambdaHistName,
                                      weakDecayFraction,
                                      lambdaXScale,
                                      sourceAlpha,
                                      shortLivedResonanceFraction,
                                      resonanceTau,
                                      FlowMode::V2YieldOnly,
                                      v2,
                                      sourceAnisotropy,
                                      primaryResonanceFile,
                                      resonanceResonanceFile);

  TH1D* hV2Geometry = BuildResonanceCurve("hCATS_v2YieldPlusGeometry",
                                          identical,
                                          withCoulomb,
                                          sourceSize,
                                          numMomBins,
                                          kMin,
                                          kMax,
                                          lambdaFileName,
                                          lambdaHistName,
                                          weakDecayFraction,
                                          lambdaXScale,
                                          sourceAlpha,
                                          shortLivedResonanceFraction,
                                          resonanceTau,
                                          FlowMode::V2YieldAndGeometry,
                                          v2,
                                          sourceAnisotropy,
                                          primaryResonanceFile,
                                          resonanceResonanceFile);

  if(!hNoFlow || !hV2Only || !hV2Geometry) {
    std::cout << "ERROR: failed to build one of the flow comparison curves.\n";
    delete hNoFlow;
    delete hV2Only;
    delete hV2Geometry;
    return;
  }

  hNoFlow->SetLineColor(kBlack);
  hNoFlow->SetLineWidth(2);
  hV2Only->SetLineColor(kBlue + 1);
  hV2Only->SetLineWidth(2);
  hV2Only->SetLineStyle(2);
  hV2Geometry->SetLineColor(kRed + 1);
  hV2Geometry->SetLineWidth(2);

  TCanvas canvas("cFlowCompare", "Flow modulation comparison", 900, 650);
  hNoFlow->Draw("hist l");
  hV2Only->Draw("hist l same");
  hV2Geometry->Draw("hist l same");

  TLegend legend(0.50, 0.68, 0.88, 0.88);
  legend.SetBorderSize(0);
  legend.SetFillStyle(0);
  legend.AddEntry(hNoFlow, "no flow", "l");
  legend.AddEntry(hV2Only, Form("v_{2}=%.3f yield only", v2), "l");
  legend.AddEntry(hV2Geometry,
                  Form("v_{2}=%.3f + s_{2}=%.2f geometry", v2, sourceAnisotropy),
                  "l");
  legend.Draw();

  TFile outputFile(outputFileName, "RECREATE");
  hNoFlow->Write("", TObject::kOverwrite);
  hV2Only->Write("", TObject::kOverwrite);
  hV2Geometry->Write("", TObject::kOverwrite);
  canvas.Write("", TObject::kOverwrite);
  outputFile.Close();

  TString pdfName = outputFileName;
  TString pngName = outputFileName;
  canvas.SaveAs(pdfName.ReplaceAll(".root", ".pdf"));
  canvas.SaveAs(pngName.ReplaceAll(".root", ".png"));

  std::cout << "Saved flow comparison to " << outputFileName << "\n";

  delete hNoFlow;
  delete hV2Only;
  delete hV2Geometry;
}
