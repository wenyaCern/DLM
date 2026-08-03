#ifndef FITCATS3D_RUNNER_CONFIG_H
#define FITCATS3D_RUNNER_CONFIG_H

#include "FitCATS3D.h"

inline FitCATS3DConfig MakeFitCATS3DBaseConfig()
{
    FitCATS3DConfig cfg;

    cfg.InputFileName = "../File/out_3DCF_mT_cent_qn_phi_tofOnly.root";
    cfg.SameHistoName = "H3_same_cent1_mt4_qn0_phi0";
    cfg.MixedHistoName = "H3_mixed_cent1_mt4_phi0";
    cfg.OutputFileName = "../fit_CATS3D_gaussian_cent1_mt4_wtCoulomb.root";

    cfg.UseSameMixed = true;
    cfg.Identical = true;
    cfg.WithCoulomb = true;
    cfg.UseLambdaFeedDown = true;
    cfg.LambdaFileName = "../File/lam/kStarDepLambda_cent10to30_pion_PbPb.root";
    cfg.LambdaHistoName = "hPurityKstar_mT0.60to0.70_cent10to30";
    cfg.weakDecayFraction = 0.12;
    cfg.RebinFactor = 1;
    cfg.MinuitMaxCalls = 20000;
    cfg.RunMinuitSimplexFirst = true;
    cfg.RunMinuitHesse = true;
    cfg.MinuitSimplexMaxCalls = 3000;
    cfg.MinuitTolerance = 0.1;
    cfg.MinuitRadiusStep = 0.4;
    cfg.UseModelCache = false;
    cfg.DebugQSymmetry = true;
    cfg.UseFastKernel = true;
    cfg.FastKernelSubdivisions = 1;
    cfg.UseExternalFastKernel = true;
    cfg.ExternalFastKernelByQ = true;
    cfg.ExternalFastKernelFileName = "../OutputFiles/dataKernel_qData_r31_coulomb_byQ.root";
    cfg.ExternalFastKernelTreeName = "CATS3DFastKernel6DByQ";
    cfg.MaxModelCacheEntries = 0;
    cfg.ModelCacheToleranceFm = 0.0;

    cfg.fitQmaxMeV = 100.0;
    cfg.normQminMeV = 100.0;
    cfg.normQmaxMeV = 200.0;
    cfg.NormalizeDataInNormWindow = true;
    cfg.NormalizeModelInNormWindow = true;
    cfg.FixNorm = true;
    cfg.UseBackground = true;
    cfg.PrefitBackground = true;
    cfg.FixBackground = true;
    cfg.BackgroundNonFemtoQminMeV = 100.0;
    cfg.sliceCutAbs = 0.04;
    cfg.relativeMomentumFactor = 0.5;
    cfg.NumGridBins = 31;

    cfg.UseResonanceSource = false;
    cfg.ResonanceWeight1 = 0.0;
    cfg.ResonanceWeight2 = 0.0;
    cfg.ResonancePRFileName = "../EPOS_data/resonance_PR.root";
    cfg.ResonancePRTreeName = "EPOS_PrimReson_LCMS";
    cfg.ResonanceRRFileName = "../EPOS_data/resonance_RR.root";
    cfg.ResonanceRRTreeName = "EPOS_Reson_LCMS";
    cfg.ResonanceDefaultTauFm = 1.5;
    cfg.ResonanceMcIter = 50000;
    cfg.ResonanceRadBins = 51;
    cfg.ResonanceAlphaInit = 2.0;
    cfg.ResonanceAlphaMin = 0.8;
    cfg.ResonanceAlphaMax = 10.0;
    cfg.ResonanceAlphaBins = 9;
    cfg.ResonanceUseEllipsoidalCore = true;
    cfg.FixResonanceAlpha = true;

    // cent1
    // cfg.RoutInit = 5.61466;
    // cfg.RsideInit = 4.90296;
    // cfg.RlongInit = 4.78343;
    cfg.RoutInit = 4.0;
    cfg.RsideInit = 4.0;
    cfg.RlongInit = 3.6;
    cfg.NormInit = 1.0;
    return cfg;
}

inline FitCATS3DConfig MakeFitCATS3DResonanceConfig()
{
    FitCATS3DConfig cfg = MakeFitCATS3DBaseConfig();
    cfg.OutputFileName = "../fit_CATS3D_resonance_cent1_mt4_wtCoulomb_levyAlpha2.0.root";
    cfg.UseResonanceSource = true;
    cfg.ResonanceWeight1 = 0.72;
    cfg.ResonanceWeight2 = 0.72;
    cfg.ResonanceAlphaInit = 2.0;
    cfg.ResonanceUseEllipsoidalCore = true;
    cfg.FixResonanceAlpha = true;
    return cfg;
}

#endif
