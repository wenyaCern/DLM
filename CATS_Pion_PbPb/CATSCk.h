#ifndef CATSCK_H
#define CATSCK_H

#include "TString.h"

void DrawCATSCk(bool identical=true,
                bool withCoulomb=true,
                double sourceSize=12.0,
                unsigned numMomBins=500,
                double kMin=0.0,
                double kMax=400.0,
                TString outputFileName="./out/CATS_ck_only.root",
                TString histName="hCATS_Ck_pipi",
                TString lambdaFileName="../Files/lam/kStarDepLambda_cent0to10_pion_PbPb.root",
                TString lambdaHistName="hLambdaKstar_mT0.30to0.50_cent0to10",
                double weakDecayFraction=0.12,
                double lambdaXScale=1000.0,
                bool useResonanceSource=false,
                double sourceAlpha=2.0,
                double shortLivedResonanceFraction=0.72,
                double resonanceTau=1.5,
                TString primaryResonanceFile="../Files/ForMax_pi_piReso_noOmega.root",
                TString resonanceResonanceFile="../Files/ForMax_piReso_piReso_noOmega.root");

void DrawCATSCkDoubleGaussianReso(bool identical=true,
                                  bool withCoulomb=true,
                                  double sourceCoreSize=12.0,
                                  double sourceTailSize=45.0,
                                  double sourceCoreFraction=0.85,
                                  unsigned numMomBins=500,
                                  double kMin=0.0,
                                  double kMax=400.0,
                                  TString outputFileName="./out/CATS_ck_double_gauss_reso.root",
                                  TString histName="hCATS_Ck_pipi_double_gauss",
                                  TString lambdaFileName="../Files/lam/kStarDepLambda_cent0to10_pion_PbPb.root",
                                  TString lambdaHistName="hLambdaKstar_mT0.30to0.50_cent0to10",
                                  double weakDecayFraction=0.12,
                                  double lambdaXScale=1000.0,
                                  double shortLivedResonanceFraction=0.72,
                                  double resonanceTau=1.5,
                                  TString primaryResonanceFile="../Files/ForMax_pi_piReso_noOmega.root",
                                  TString resonanceResonanceFile="../Files/ForMax_piReso_piReso_noOmega.root");

#endif
