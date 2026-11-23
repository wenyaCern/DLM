#include "FitCATS3DRunnerConfig.h"

#include "TROOT.h"

#include <filesystem>
#include <iostream>
#include <stdexcept>

int main(int argc, char** argv)
{
    gROOT->SetBatch(true);
    if(argc < 4){
        std::cerr << "Usage: fit_CATS3D_derived_full MT_INDEX QN_INDEX PHI_INDEX [axis-slices] [rinit=R] [rmin=R] [fitqmax=MEV] [lambda-after-dca-cut] [weak155] [adaptive-kernel=PATH] [smearfile=PATH] [output=PATH]\n";
        return 2;
    }
    try{
        auto cfg = MakeFitCATS3DDerivedFullBminusConfig(
            std::stoi(argv[1]), std::stoi(argv[2]), std::stoi(argv[3]));
        for(int i=4; i<argc; ++i){
            const std::string arg(argv[i]);
            if(arg == "axis-slices"){
                cfg.FitAxisSlicesOnly = true;
                cfg.OutputFileName.ReplaceAll("_outerbox_weak14_fit/", "_outerbox_weak14_axis_slices_fit/");
            }
            else if(arg == "lambda083"){
                cfg.LambdaFileName = "derivedFull_Bminus_cent0to10_anyCPR_lambda083_test/lambda083.root";
                cfg.LambdaHistoName = "lambda083";
                cfg.LambdaMtBinIndex = -1;
                cfg.OutputFileName.ReplaceAll("_outerbox_weak14_fit/", "_lambda083_weak14_fit/");
                cfg.OutputFileName.ReplaceAll(".root", "_lambda083.root");
            }
            else if(arg == "lambda085"){
                cfg.LambdaFileName = "derivedFull_Bminus_cent0to10_anyCPR_lambda083_test/lambda085.root";
                cfg.LambdaHistoName = "lambda085";
                cfg.LambdaMtBinIndex = -1;
                cfg.OutputFileName.ReplaceAll(".root", "_lambda085.root");
            }
            else if(arg == "fine-kernel-derived-lambda"){
                cfg.LambdaFileName = CATS3DPathFromEnv(
                    "CATS3D_LAMBDA_FILE",
                    "../../inputs/data/PbPb_cent0to10_lambda_gen_derived_fullPt.root");
                cfg.LambdaHistoName = "lambdaFitWindow";
                cfg.UseFineModelGrid = true;
                cfg.ExternalFastKernelFileName = CATS3DPathFromEnv(
                    "CATS3D_KERNEL_FILE", "../../inputs/OutputFiles/sourceKernel_q10_r31_rmax20_coulomb_byQ.root");
                cfg.gridMinFm = -20.0;
                cfg.gridMaxFm = 20.0;
                cfg.BackgroundNonFemtoQminMeV = 100.0;
                cfg.BackgroundNonFemtoQmaxMeV = 200.0;
                cfg.OutputFileName.ReplaceAll("_outerbox_weak14_fit/", "_fineKernel_derivedLambda_baseline100to200_fit/");
            }
            else if(arg.rfind("adaptive-kernel=",0) == 0){
                cfg.UseAdaptiveKernel=true;
                cfg.AdaptiveKernelFileName=arg.substr(16).c_str();
                cfg.UseExternalFastKernel=false;
                cfg.UseFineModelGrid=false;
                cfg.UseModelCache=true;
                cfg.ModelCacheToleranceFm=0.02;
                cfg.OutputFileName.ReplaceAll("_outerbox_weak14_fit/", "_adaptiveKernel_outerbox_weak14_fit/");
            }
            else if(arg == "derived-lambda"){
                cfg.LambdaFileName = CATS3DPathFromEnv(
                    "CATS3D_LAMBDA_FILE",
                    "../../inputs/data/PbPb_cent0to10_lambda_gen_derived_fullPt.root");
                cfg.LambdaHistoName = "lambdaFitWindow";
            }
            else if(arg == "lambda-after-dca-cut"){
                cfg.LambdaHistoName = "lambdaAfterDCACut";
            }
            else if(arg == "weak155"){
                cfg.weakDecayFraction = 0.155;
            }
            else if(arg.rfind("weakfraction=",0) == 0){
                const double fraction = std::stod(arg.substr(13));
                if(fraction < 0.0 || fraction >= 1.0)
                    throw std::runtime_error("weakfraction must be within [0, 1)");
                cfg.weakDecayFraction = fraction;
            }
            else if(arg == "smear3d"){
                cfg.UseFactorizedQSmearing = true;
                cfg.QSmearingFileName = CATS3DPathFromEnv(
                    "CATS3D_SMEARING_FILE",
                    "../../inputs/data/AnalysisResults_lambdaPara_smearing3DLCMS_test.root");
            }
            else if(arg.rfind("smearfile=",0) == 0){
                cfg.UseFactorizedQSmearing = true;
                cfg.QSmearingFileName = arg.substr(10).c_str();
            }
            else if(arg == "diagonal-band"){
                cfg.UseDiagonalBandBackground = true;
                // B- pi+pi+: measured non-femto stripe follows
                // q_side ~= 0.15 q_out in this low-mT bin.
                cfg.BandSlopeSidePerOut = 0.15;
                cfg.BandWidthPerpGeV = 0.030;
                cfg.BandWidthLongGeV = 0.030;
                cfg.BandAmplitudeInit = 0.0103;
            }
            else if(arg == "band-fixed"){
                cfg.FixBandAmplitude = true;
            }
            else if(arg == "kernel-fitbox-only"){
                cfg.OutputBaselineOutsideFitBox = true;
            }
            else if(arg.rfind("output=",0) == 0){
                cfg.OutputFileName = arg.substr(7).c_str();
            }
            else if(arg.rfind("rinit=",0) == 0){
                const double r = std::stod(arg.substr(6));
                if(r < cfg.RadiusMinFm || r > cfg.RadiusMaxFm)
                    throw std::runtime_error("rinit outside radius bounds");
                cfg.RoutInit = cfg.RsideInit = cfg.RlongInit = r;
                TString tag = Form("_rinit%.1f.root", r);
                tag.ReplaceAll(".", "p");
                tag.ReplaceAll("proot", ".root");
                cfg.OutputFileName.ReplaceAll(".root", tag);
            }
            else if(arg.rfind("rmin=",0) == 0){
                const double r = std::stod(arg.substr(5));
                if(r < cfg.RadiusMinFm || r >= cfg.RadiusMaxFm)
                    throw std::runtime_error("rmin outside radius bounds");
                cfg.RadiusMinFm = r;
                TString tag = Form("_rmin%.1f.root", r);
                tag.ReplaceAll(".", "p");
                tag.ReplaceAll("proot", ".root");
                cfg.OutputFileName.ReplaceAll(".root", tag);
            }
            else if(arg.rfind("fitqmax=",0) == 0){
                const double qmax = std::stod(arg.substr(8));
                if(qmax <= 0.0 || qmax > 200.0)
                    throw std::runtime_error("fitqmax must be within (0, 200] MeV/c");
                cfg.fitQmaxMeV = qmax;
            }
            else throw std::runtime_error("Unknown fit mode");
        }
        std::filesystem::create_directories(
            std::filesystem::path(cfg.OutputFileName.Data()).parent_path());
        const auto result = FitDataCATS3D(cfg);
        return result.FitStatus == 0 ? 0 : 3;
    }
    catch(const std::exception& e){
        std::cerr << "fit_CATS3D_derived_full: " << e.what() << '\n';
        return 1;
    }
}
