#ifndef FITCATS3D_RUNNER_CONFIG_H
#define FITCATS3D_RUNNER_CONFIG_H

#include <cstdlib>
#include <stdexcept>

#include "FitCATS3D.h"

inline TString CATS3DPathFromEnv(const char* name, const char* fallback)
{
    const char* value = std::getenv(name);
    return value && *value ? TString(value) : TString(fallback);
}

inline FitCATS3DConfig MakeFitCATS3DBaseConfig()
{
    FitCATS3DConfig cfg;

    cfg.InputFileName = "../../inputs/File/out_3DCF_mT_cent_qn_phi_tofOnly_raw.root";
    cfg.SameHistoName = "H3_same_cent1_mt7_qn0_phi0";
    cfg.MixedHistoName = "H3_mixed_cent1_mt7_phi0";
    cfg.OutputFileName = "../fit_CATS3D_gaussian_cent1_mt7_wtCoulomb.root";

    cfg.UseSameMixed = true;
    cfg.Identical = true;
    cfg.WithCoulomb = true;
    cfg.UseLambdaFeedDown = true;
    cfg.LambdaFileName = "../../inputs/File/lam/kStarDepLambda_cent10to30_pion_PbPb.root";
    cfg.LambdaHistoName = "hPurityKstar_mT1.10to1.51_cent10to30";
    // Pb-Pb weak-decay fraction estimated with Thermal-FIST.
    cfg.weakDecayFraction = 0.14;
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
    cfg.ExternalFastKernelFileName = "../../inputs/OutputFiles/dataKernel_qData_r31_coulomb_byQ.root";
    cfg.ExternalFastKernelTreeName = "CATS3DFastKernel6DByQ";
    cfg.MaxModelCacheEntries = 0;
    cfg.ModelCacheToleranceFm = 0.0;

    cfg.fitQmaxMeV = 100.0;
    // Raw-count exposure normalization in the outer Cartesian region:
    // |q_i| <= 200 MeV/c, excluding the central |q_i| < 150 MeV/c cube.
    cfg.normQminMeV = 150.0;
    cfg.normQmaxMeV = 200.0;
    cfg.NormalizeDataInNormWindow = true;
    cfg.NormalizeModelInNormWindow = true;
    cfg.UsePmlChi2 = true;
    cfg.FixNorm = true;
    cfg.UseBackground = true;
    cfg.PrefitBackground = true;
    cfg.FixBackground = true;
    cfg.BackgroundNonFemtoQminMeV = 150.0;
    cfg.sliceCutAbs = 0.04;
    cfg.relativeMomentumFactor = 0.5;
    cfg.NumGridBins = 31;

    cfg.UseResonanceSource = false;
    cfg.ResonanceWeight1 = 0.0;
    cfg.ResonanceWeight2 = 0.0;
    cfg.ResonancePRFileName = "../../inputs/EPOS_data/resonance_PR.root";
    cfg.ResonancePRTreeName = "EPOS_PrimReson_LCMS";
    cfg.ResonanceRRFileName = "../../inputs/EPOS_data/resonance_RR.root";
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
    cfg.RoutInit = 2.80021;
    cfg.RsideInit = 3.07897;
    cfg.RlongInit = 2.65327;
    cfg.NormInit = 1.0;
    return cfg;
}

inline FitCATS3DConfig MakeFitCATS3DDerivedFullBminusConfig(
    const int iMt, const int iQn, const int iPhi)
{
    if(iMt < 0 || iMt >= 3 || iQn < 0 || iQn >= 5 ||
       iPhi < 0 || iPhi >= 6){
        throw std::runtime_error("Derived full Bminus indices must be mT 0..2, qn 0..4, phi 0..5");
    }
    FitCATS3DConfig cfg = MakeFitCATS3DBaseConfig();
    const TString stem = Form("mt%d_qn%d_phi%d", iMt, iQn, iPhi);
    cfg.InputFileName = CATS3DPathFromEnv(
        "CATS3D_DERIVED_CF",
        "../../inputs/data/CF3D_derivedFull_Bminus_cent0to10_anyCPR.root");
    cfg.SameHistoName = stem + "/same3d";
    cfg.MixedHistoName = stem + "/mixed3d";
    const TString outputDir = CATS3DPathFromEnv(
        "CATS3D_FIT_OUTPUT_DIR",
        "../../fit/results/derivedFull_Bminus_cent0to10_anyCPR_outerbox_weak14_fit");
    cfg.OutputFileName = Form("%s/mt%d/qn%d/fit_phi%d.root",
                              outputDir.Data(), iMt, iQn, iPhi);

    cfg.LambdaFileName = CATS3DPathFromEnv(
        "CATS3D_LAMBDA_FILE",
        "../../inputs/data/PbPb_cent0to10_piPlus_lambda_gen_fullPt_requested_mt.root");
    cfg.LambdaHistoName = "lambdaFitWindow";
    cfg.LambdaMtBinIndex = iMt + 1;
    // Omega weak-decay correction is not included in this lambda_gen input.
    cfg.weakDecayFraction = 0.14;
    cfg.UseLambdaFeedDown = true;

    cfg.RebinFactor = 1;
    cfg.fitQmaxMeV = 100.0;
    cfg.normQminMeV = 120.0;
    cfg.normQmaxMeV = 180.0;
    cfg.UseWedgeNormalizationRegion = true;
    cfg.NormWedgeSideOffsetMeV = 40.0;
    cfg.NormalizeDataInNormWindow = true;
    cfg.NormalizeModelInNormWindow = true;

    // Fit the baseline only in the outer Cartesian box; never use the peak.
    cfg.UseBackground = true;
    cfg.PrefitBackground = true;
    cfg.UsePmlBackgroundPrefit = true;
    cfg.UseDiagonalQuadraticBackground = true;
    cfg.UseHybridBackgroundRegion = false;
    cfg.UseSphericalBackgroundRegion = false;
    cfg.BackgroundNonFemtoQminMeV = 120.0;
    cfg.BackgroundNonFemtoQmaxMeV = 200.0;
    cfg.FixBackground = true;
    cfg.FixNorm = true;
    cfg.UsePmlChi2 = true;
    cfg.sliceCutAbs = 0.02;

    cfg.MinuitMaxCalls = 8000;
    cfg.RunMinuitSimplexFirst = true;
    cfg.RunMinuitHesse = true;
    cfg.MinuitSimplexMaxCalls = 1200;
    cfg.MinuitTolerance = 0.1;
    cfg.MinuitRadiusStep = 0.3;
    const double r0[] = {6.5, 5.5, 4.5};
    cfg.RoutInit = r0[iMt];
    cfg.RsideInit = r0[iMt];
    cfg.RlongInit = r0[iMt];
    cfg.RadiusMinFm = 0.5;
    cfg.RadiusMaxFm = 16.0;
    cfg.UseResonanceSource = false;
    cfg.ResonanceAlphaInit = 2.0;
    cfg.FixResonanceAlpha = true;

    cfg.UseFastKernel = false;
    cfg.UseExternalFastKernel = false;
    cfg.UseAdaptiveKernel = true;
    cfg.AdaptiveKernelFileName = CATS3DPathFromEnv(
        "CATS3D_KERNEL_FILE", "../../inputs/OutputFiles/adaptive_pipi_coulomb_q19_r34_byQ.root");
    cfg.UseModelCache = true;
    cfg.ModelCacheToleranceFm = 0.02;
    cfg.NumGridBins = 31;
    cfg.relativeMomentumFactor = 0.5;
    return cfg;
}

inline const char* FitCATS3DLambdaMtTag(const int iMt)
{
    switch(iMt){
        case 2: return "0.40to0.50";
        case 3: return "0.50to0.60";
        case 4: return "0.60to0.70";
        case 5: return "0.70to0.90";
        case 6: return "0.90to1.10";
        case 7: return "1.10to1.51";
        default: return nullptr;
    }
}

inline FitCATS3DConfig MakeFitCATS3DMtConfig(const int iMt)
{
    FitCATS3DConfig cfg = MakeFitCATS3DBaseConfig();
    const char* mtTag = FitCATS3DLambdaMtTag(iMt);
    if(!mtTag) throw std::runtime_error("Unsupported mT index");
    cfg.SameHistoName = Form("H3_same_cent1_mt%d_qn0_phi0", iMt);
    cfg.MixedHistoName = Form("H3_mixed_cent1_mt%d_phi0", iMt);
    cfg.LambdaHistoName = Form("hPurityKstar_mT%s_cent10to30", mtTag);
    cfg.OutputFileName = Form(
        "../fit_CATS3D_gaussian_cent1_mt%d_wtCoulomb.root", iMt);
    return cfg;
}

inline FitCATS3DConfig MakeFitCATS3DSubSampleConfig(const int iSamp)
{
    FitCATS3DConfig cfg = MakeFitCATS3DBaseConfig();
    cfg.InputFileName = Form(
        "../../inputs/File/subSamp/out_3DCF_mT_cent_qn_phi_tofOnly_samp%d.root", iSamp);
    cfg.OutputFileName = Form(
        "../fit_CATS3D_gaussian_cent1_mt7_samp%d.root", iSamp);
    cfg.LoadBackgroundShape = true;
    cfg.BackgroundShapeFileName =
        "../fit_CATS3D_gaussian_cent1_mt7_wtCoulomb.root";
    cfg.PrefitBackground = true;
    cfg.FixBackground = true;
    cfg.FixNorm = false;
    return cfg;
}

inline FitCATS3DConfig MakeFitCATS3DBootstrapConfig(const int iReplica)
{
    // A bootstrap replica must repeat the same estimator used for the full
    // data.  In particular, prefit all ten baseline parameters on this
    // replica, then keep that replica-specific baseline fixed during the
    // source fit.  Loading the full-data baseline would condition away its
    // statistical fluctuations and underestimate the radius uncertainty.
    FitCATS3DConfig cfg = MakeFitCATS3DBaseConfig();
    cfg.InputFileName = Form(
        "../../inputs/File/bootstrap/out_3DCF_mT_cent_qn_phi_tofOnly_boot%d.root", iReplica);
    cfg.OutputFileName = Form(
        "../fit_CATS3D_gaussian_cent1_mt7_boot_fullBaseline%d.root", iReplica);
    cfg.LoadBackgroundShape = false;
    cfg.BackgroundShapeFileName = "";
    cfg.PrefitBackground = true;
    cfg.FixBackground = true;
    cfg.FixNorm = true;
    return cfg;
}

inline FitCATS3DConfig MakeFitCATS3DBootstrapMtConfig(const int iMt,
                                                      const int iReplica)
{
    FitCATS3DConfig cfg = MakeFitCATS3DMtConfig(iMt);
    cfg.InputFileName = Form(
        "../../inputs/File/bootstrap/mt%d/out_3DCF_mT_cent_qn_phi_tofOnly_boot%d.root",
        iMt, iReplica);
    cfg.OutputFileName = Form(
        "../fit_CATS3D_gaussian_cent1_mt%d_boot_fullBaseline%d.root",
        iMt, iReplica);
    cfg.LoadBackgroundShape = false;
    cfg.BackgroundShapeFileName = "";
    cfg.PrefitBackground = true;
    cfg.FixBackground = true;
    cfg.FixNorm = true;
    return cfg;
}

inline FitCATS3DConfig MakeFitCATS3DNominalBootstrapMtConfig(const int iMt,
                                                             const int iReplica)
{
    FitCATS3DConfig cfg = MakeFitCATS3DMtConfig(iMt);
    cfg.InputFileName = Form(
        "../../inputs/File/bootstrap/mt%d/out_3DCF_mT_cent_qn_phi_tofOnly_boot%d.root",
        iMt, iReplica);
    cfg.OutputFileName = Form(
        "../fit_CATS3D_gaussian_cent1_mt%d_boot_nominalHybridOldLambdaWeak12_%d.root",
        iMt, iReplica);
    cfg.UseDiagonalQuadraticBackground = true;
    cfg.UseHybridBackgroundRegion = true;
    cfg.UsePmlBackgroundPrefit = true;
    cfg.BackgroundNonFemtoQminMeV = 100.0;
    cfg.BackgroundNonFemtoQmaxMeV = 200.0;
    cfg.sliceCutAbs = 0.04;
    cfg.LoadBackgroundShape = false;
    cfg.PrefitBackground = true;
    cfg.FixBackground = true;
    cfg.FixNorm = true;
    cfg.LambdaFileName =
        "../../inputs/File/lam/lambda_LHC24g3/kStarDepLambda_cent10to30_pion_PbPb.root";
    cfg.weakDecayFraction = 0.12;
    return cfg;
}

inline FitCATS3DConfig MakeFitCATS3DResonanceConfig()
{
    FitCATS3DConfig cfg = MakeFitCATS3DBaseConfig();
    cfg.OutputFileName = "../fit_CATS3D_resonance_cent1_mt7_wtCoulomb_levyAlpha2.0.root";
    cfg.UseResonanceSource = true;
    cfg.ResonanceWeight1 = 0.72;
    cfg.ResonanceWeight2 = 0.72;
    cfg.ResonanceAlphaInit = 2.0;
    cfg.ResonanceUseEllipsoidalCore = true;
    cfg.FixResonanceAlpha = true;
    return cfg;
}

#endif
