#ifndef FITCATS3D_H
#define FITCATS3D_H

#include "TString.h"
#include <memory>

class TH3D;
class DLM_CleverReso3D;

struct FitCATS3DConfig {
    TString InputFileName = "./File/out_3DCF_mT_cent_qn_phi_tofOnly.root";
    TString SameHistoName = "H3_same_cent0_mt3_qn0_phi0";
    TString MixedHistoName = "H3_mixed_cent0_mt3_phi0";
    TString CfHistoName = "";
    TString OutputFileName = "./fit_CATS3D_gaussian.root";
    TString LambdaFileName = "../CATS_Pion_PbPb/Files/lam/kStarDepPurity_cent0to10_pion_PbPb.root";
    TString LambdaHistoName = "hPurityKstar_mT0.50to0.60_cent0to10";
    // Select a physical mT row when LambdaHistoName points to a TH2.
    int LambdaMtBinIndex = 0;

    bool UseSameMixed = true;
    bool Identical = true;
    bool WithCoulomb = true;
    bool UseLambdaFeedDown = false;
    bool UseResonanceSource = false;

    int RebinFactor = 2;
    int MinuitMaxCalls = 80;
    bool RunMinuitSimplexFirst = false;
    bool RunMinuitHesse = false;
    int MinuitSimplexMaxCalls = 2000;
    double MinuitTolerance = 1.0;
    double MinuitRadiusStep = 0.1;
    bool UseModelCache = true;
    bool EnableProfiling = true;
    bool DebugQSymmetry = false;
    bool UseFastKernel = true;
    unsigned FastKernelSubdivisions = 1;
    // Evaluate CATS on a 10-MeV grid inside the fit box, then average each
    // eight fine cells into the corresponding 20-MeV data cell.
    bool UseFineModelGrid = false;
    bool UseExternalFastKernel = false;
    // Nonuniform positive-q wavefunction kernel, folded with the Gaussian source.
    bool UseAdaptiveKernel = false;
    TString AdaptiveKernelFileName = "";
    bool UseFactorizedQSmearing = false;
    TString QSmearingFileName = "";
    TString QSmearingHistogramBase = "femto-pair-track-track/TrackTrack/SE/MC/";
    bool ExternalFastKernelByQ = true;
    TString ExternalFastKernelFileName = "";
    TString ExternalFastKernelTreeName = "CATS3DFastKernel6DByQ";
    unsigned MaxModelCacheEntries = 12;
    double ModelCacheToleranceFm = 0.08;

    double qAxisScaleToMeV = 1000.0;
    double relativeMomentumFactor = 0.5;
    double lambdaAxisScaleToMeV = 1000.0;
    double weakDecayFraction = 0.12;
    double fitQmaxMeV = 200.0;
    bool FitAxisSlicesOnly = false;
    double FitAxisSliceCutMeV = 20.0;
    double normQminMeV = 300.0;
    double normQmaxMeV = 500.0;
    bool UseSphericalNormalizationRegion = false;
    // Spherical shell in |q| with the out-side diagonal band excluded.
    bool UseWedgeNormalizationRegion = false;
    double NormWedgeSideOffsetMeV = 40.0;
    bool UseRadialQuadraticBackground = false;
    bool UseDiagonalQuadraticBackground = false;
    // Residual close-pair depletion aligned with the qside-qout band.
    bool UseDiagonalBandBackground = false;
    double BandSlopeSidePerOut = 0.0;
    double BandWidthPerpGeV = 0.03;
    double BandWidthLongGeV = 0.03;
    double BandAmplitudeInit = 0.01;
    bool FixBandAmplitude = false;
    // In saved whole-range maps, display the prefit baseline alone outside
    // the source-fit box; the kernel is only constrained inside that box.
    bool OutputBaselineOutsideFitBox = false;
    bool UseDiagonalQuarticBackground = false;
    bool UseSphericalBackgroundRegion = false;
    bool UseHybridBackgroundRegion = false;
    bool FloatDiagonalQuadraticBackgroundInFinalFit = false;
    bool NormalizeDataInNormWindow = true;
    bool NormalizeModelInNormWindow = true;
    bool UsePmlChi2 = true;
    bool FixNorm = true;
    bool UseBackground = true;
    bool PrefitBackground = true;
    bool UsePmlBackgroundPrefit = false;
    bool FixBackground = true;
    bool LoadBackgroundShape = false;
    TString BackgroundShapeFileName = "";
    double BackgroundNonFemtoQminMeV = 100.0;
    double BackgroundNonFemtoQmaxMeV = 200.0;
    double sliceCutAbs = 0.06;

    unsigned NumGridBins = 21;
    double gridMinFm = -50.0;
    double gridMaxFm = 50.0;

    TString ResonancePRFileName = "./EPOS_data/resonance_PR.root";
    TString ResonancePRTreeName = "EPOS_PrimReson_LCMS";
    TString ResonanceRRFileName = "./EPOS_data/resonance_RR.root";
    TString ResonanceRRTreeName = "EPOS_Reson_LCMS";
    bool UseResonancePR = true;
    bool UseResonanceRR = true;
    double ResonanceWeight1 = 0.0;
    double ResonanceWeight2 = 0.0;
    double ResonanceDefaultTauFm = 1.5;
    unsigned ResonanceMaxPairs = 0;
    unsigned ResonanceMcIter = 200000;
    unsigned ResonanceRadBins = 81;
    double ResonanceRadMinFm = -64.0;
    double ResonanceRadMaxFm = 64.0;
    unsigned ResonanceScaleBins = 128;
    double ResonanceScaleMinFm = 0.2;
    double ResonanceScaleMaxFm = 30.0;
    double ResonanceAlphaInit = 1.4;
    double ResonanceAlphaMin = 1.2;
    double ResonanceAlphaMax = 1.6;
    unsigned ResonanceAlphaBins = 9;
    bool FixResonanceAlpha = false;
    bool ResonanceUseEllipsoidalCore = true;
    bool ResonanceDebug = false;

    double RoutInit = 9.0;
    double RsideInit = 9.0;
    double RlongInit = 9.0;
    double RadiusMinFm = 1.0;
    double RadiusMaxFm = 16.0;
    double NormInit = 1.0;
};

struct FitCATS3DResult {
    double Rout = 0.0;
    double Rside = 0.0;
    double Rlong = 0.0;
    double Alpha = 2.0;
    double RoutErr = 0.0;
    double RsideErr = 0.0;
    double RlongErr = 0.0;
    double AlphaErr = 0.0;
    double Norm = 1.0;
    double BandAmplitude = 0.0;
    double BandAmplitudeErr = 0.0;
    double Chi2 = 0.0;
    int Ndf = 0;
    int NFitBins = 0;
    int FitStatus = -1;
    int MinuitIstat = -1;
    int FitQuality = 3;
    double MinuitAmin = 0.0;
    double MinuitEdm = 0.0;
    double MinuitErrDef = 1.0;
    int MinuitNvpar = 0;
    int MinuitNparx = 0;
    double BkgPar[10] = {1.0, 0.0, 0.0, 0.0, 0.0,
                         0.0, 0.0, 0.0, 0.0, 0.0};
    double BkgParErr[10] = {0.0, 0.0, 0.0, 0.0, 0.0,
                            0.0, 0.0, 0.0, 0.0, 0.0};

    double Chi2Ndf() const { return Ndf > 0 ? Chi2 / double(Ndf) : 0.0; }
};

TH3D* BuildDataCk3D(const FitCATS3DConfig& cfg);
FitCATS3DResult FitDataCATS3D(const FitCATS3DConfig& cfg = FitCATS3DConfig());
FitCATS3DResult FitDataCATS3D(const FitCATS3DConfig& cfg,
                              DLM_CleverReso3D* resonanceSource);
std::unique_ptr<DLM_CleverReso3D> BuildCATS3DResonanceSource(const FitCATS3DConfig& cfg);
void fitDataCATS3D();

#endif
