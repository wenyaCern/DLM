#ifndef BASICS3D_H
#define BASICS3D_H

#include "/Users/wuwenya/Documents/code/femto_flow/CATS/CATS_install/install/include/CATS3D.h"
#include "TH3D.h"
#include "TString.h"

double Basics_Source_Gauss3D(double* Pars);

void Basics_PiPiCATS3D_Setup(CATS3D& Kitty,
                             const bool& Identical=true,
                             const bool& WithCoulomb=true);

TH3D* Basics_CATS3D_GetCkHisto(CATS3D& Kitty,
                               const TString& HistoName="hCk3D");

void Save_TH3D_to_ROOT(const TString& FileName,
                       TH3D* histo,
                       const TString& HistoName="");

void Basics_RunSourceComparison3D(const TString& outFileName="OutputFiles/source_compare_R5.root");

void Basics_RunGaussianR20DataRangeCheck3D(
    const TString& outFileName="OutputFiles/gauss_R20_data_range_check.root",
    const double qAbsMaxMeV=1000.0,
    const unsigned nQBins=21,
    const bool runCoulomb=true,
    const bool runNoCoulomb=true);

void Basics_RunFineGaussianNoCoulomb3D(
    const TString& outFileName="OutputFiles/fine_gaussian_noCoulomb_R6.root",
    const double radiusFm=6.0,
    const double qAbsMaxMeV=200.0,
    const unsigned nQBins=101,
    const double sliceCutMeV=40.0,
    const bool useFastKernel=true,
    const unsigned fastSubdivisions=1);

#endif
