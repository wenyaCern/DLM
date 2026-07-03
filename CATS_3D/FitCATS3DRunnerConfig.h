#ifndef FITCATS3D_RUNNER_CONFIG_H
#define FITCATS3D_RUNNER_CONFIG_H

#include "FitCATS3D.h"

inline FitCATS3DConfig MakeFitCATS3DBaseConfig()
{
    FitCATS3DConfig cfg;

    cfg.InputFileName = "../File/out_3DCF_mT_cent_qn_phi_tofOnly.root";
    cfg.SameHistoName = "H3_same_cent0_mt2_qn0_phi0";
    cfg.MixedHistoName = "H3_mixed_cent0_mt2_phi0";
    cfg.OutputFileName = "../fit_CATS3D_gaussian_cent0_mt2_noCoulomb.root";

    cfg.UseSameMixed = true;
    cfg.Identical = true;
    cfg.WithCoulomb = true;
    cfg.UseLambdaFeedDown = true;
    cfg.LambdaFileName = "../File/lam/kStarDepLambda_cent0to10_pion_PbPb.root";
    cfg.LambdaHistoName = "hPurityKstar_mT0.40to0.50_cent0to10";
    cfg.weakDecayFraction = 0.12;
    cfg.RebinFactor = 1;
    cfg.MinuitMaxCalls = 1000;
    cfg.UseModelCache = true;
    cfg.UseFastKernel = true;
    cfg.FastKernelSubdivisions = 2;
    cfg.MaxModelCacheEntries = 12;
    cfg.ModelCacheToleranceFm = 0.08;

    cfg.fitQmaxMeV = 200.0;
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
    cfg.NumGridBins = 15;

    cfg.UseResonanceSource = false;
    cfg.ResonanceWeight1 = 0.0;
    cfg.ResonanceWeight2 = 0.0;
    cfg.ResonancePRFileName = "../EPOS_data/resonance_PR.root";
    cfg.ResonancePRTreeName = "EPOS_PrimReson_LCMS";
    cfg.ResonanceRRFileName = "../EPOS_data/resonance_RR.root";
    cfg.ResonanceRRTreeName = "EPOS_Reson_LCMS";
    cfg.ResonanceDefaultTauFm = 1.5;
    cfg.ResonanceMcIter = 20000;
    cfg.ResonanceRadBins = 51;
    cfg.ResonanceAlphaInit = 1.0;
    cfg.ResonanceAlphaMin = 0.8;
    cfg.ResonanceAlphaMax = 2.0;
    cfg.ResonanceAlphaBins = 9;
    cfg.FixResonanceAlpha = true;

    cfg.RoutInit = 8.0;
    cfg.RsideInit = 8.0;
    cfg.RlongInit = 8.0;
    cfg.NormInit = 1.0;

    return cfg;
}

inline FitCATS3DConfig MakeFitCATS3DResonanceConfig()
{
    FitCATS3DConfig cfg = MakeFitCATS3DBaseConfig();
    cfg.OutputFileName = "../fit_CATS3D_resonance_cent0_mt2_noCoulomb.root";
    cfg.UseResonanceSource = true;
    cfg.ResonanceWeight1 = 0.72;
    cfg.ResonanceWeight2 = 0.72;
    cfg.ResonanceAlphaInit = 1.6;
    cfg.FixResonanceAlpha = true;
    return cfg;
}

#endif
