#ifndef PION_PBPB_EXTENDED_CK_H
#define PION_PBPB_EXTENDED_CK_H

#include "TString.h"

class DLM_Ck;
class DLM_CkDecomposition;

extern double sourceAlpha;

void Init_pipi_Theory(TString SourceType,
                      const bool Identical,
                      const bool WithCoulomb,
                      const bool onlyTheory,
                      const double SourceSize,
                      const double SourceAlpha=2.0,
                      const unsigned NumMomBins=500,
                      const double kMin=0,
                      const double kMax=400,
                      TString lambdaFileName="../Files/lam/kStarDepLambda_cent0to10_pion_PbPb.root",
                      TString lambdaHistName="hPurityKstar_mT0.40to0.50_cent0to10");

void fitDataCoreReso(bool Identical=true,
                     bool WithCoulomb=true,
                     double alphaInit=2.0,
                     bool fixAlpha=true);

void fitDataCoreResoFixedAlpha(bool Identical=true,
                               bool WithCoulomb=true,
                               double alpha=2.0);

void fitDataCoreResoFreeAlpha(bool Identical=true,
                              bool WithCoulomb=true,
                              double alphaInit=1.5);

// Compatibility aliases kept for old main macros. They now use only
// DLM_CleverMcLevyResoTM; other source families were moved to backup.
void fitData(bool Identical=true, bool WithCoulomb=true);
void fitDataGaussOnly(bool Identical=true, bool WithCoulomb=true);
void fitData3DTo1D(bool Identical=true,
                   bool WithCoulomb=true,
                   bool usePRF3DShift=false);

void Ck_pL_Ledni_Usmani();
void RootFile_DlmCk(const TString& RootFileName,
                    const TString& GraphName,
                    DLM_Ck* CkToPlot);
void RootFile_DlmCk(const TString& RootFileName,
                    const TString& GraphName,
                    DLM_CkDecomposition* CkToPlot,
                    bool PlotIndividualContributions=false);
void Ck_pp_Decomposition(const TString& SourceType);

#endif
