#include "TString.h"

class DLM_Ck;
class DLM_CkDecomposition;

extern double sourceAlpha;

void Init_pipi_Theory(TString SourceType,
                      const bool Identical,
                      const bool WithCoulomb,
                      const bool onlyTheory,
                      const double SourceSize,
                      const double SourceAlpha=1.5,
                      const unsigned NumMomBins=500,
                      const double kMin=0,
                      const double kMax=400,
                      TString lambdaFileName="../Files/kStarDepLambda_pion_PbPb.root",
                      TString lambdaHistName="hLambdaKstar_mT0.30to0.50_cent0to10");
void fitData(bool Identical=true, bool WithCoulomb=true);
void fitData3DTo1D(bool Identical=true,
                   bool WithCoulomb=true,
                   bool usePRF3DShift=false);
void ScanAlphaBinning(bool Identical=true,
                      bool WithCoulomb=true,
                      double alphaMin=1.0,
                      double alphaMax=2.0,
                      double alphaStep=0.2,
                      double fitMin=2.0,
                      double fitMax=200.0,
                      double baselineFitMin=200.0,
                      double baselineFitMax=400.0,
                      unsigned catsMomBins=500);
void ScanAlpha8MeV(bool Identical=true,
                   bool WithCoulomb=true,
                   double alphaMin=1.0,
                   double alphaMax=2.0,
                   double alphaStep=0.2,
                   double fitMin=2.0,
                   double fitMax=200.0,
                   double baselineFitMin=200.0,
                   double baselineFitMax=400.0,
                   unsigned catsMomBins=500);
void ScanWeakDecay8MeV(bool Identical=true,
                       bool WithCoulomb=true,
                       double fixedAlpha=1.5,
                       double weakMin=0.12,
                       double weakMax=0.18,
                       double weakStep=0.02,
                       double fitMin=2.0,
                       double fitMax=200.0,
                       double baselineFitMin=200.0,
                       double baselineFitMax=400.0,
                       unsigned catsMomBins=500);

void Ck_pL_Ledni_Usmani();
//save in a file the correlation function from DLM_Ck
void RootFile_DlmCk(const TString& RootFileName, const TString& GraphName, DLM_Ck* CkToPlot);
//save in a file the correlation function from DLM_CkDecomposition
void RootFile_DlmCk(const TString& RootFileName, const TString& GraphName, DLM_CkDecomposition* CkToPlot, bool PlotIndividualContributions=false);

void Ck_pp_Decomposition(const TString& SourceType);
