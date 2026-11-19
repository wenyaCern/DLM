#include "FitCATS3D.h"
#include "FitCATS3DRunnerConfig.h"

#include "DLM_Source.h"

#include "TROOT.h"

#include <filesystem>
#include <iostream>
#include <memory>
#include <string>

int main(int argc, char** argv)
{
    gROOT->SetBatch(true);
    const bool gaussianControl = argc == 2 && std::string(argv[1]) == "gaussian-control";
    if(argc > 2 || (argc == 2 && !gaussianControl)){
        std::cerr << "Usage: fit_resonance6291_mt0 [gaussian-control]\n";
        return 2;
    }
    auto cfg = MakeFitCATS3DDerivedFullBminusConfig(0, 0, 0);
    cfg.OutputFileName = "derivedFull_Bminus_cent0to10_anyCPR_resonance6291_q10_r20_weak1557_fit/mt0/qn0/fit_phi0.root";
    cfg.LambdaFileName = "derivedFull_Bminus_cent0to10_anyCPR_lambda083_test/lambda083.root";
    cfg.LambdaHistoName = "lambda083";
    cfg.LambdaMtBinIndex = -1;
    cfg.weakDecayFraction = 0.1557;

    cfg.UseFineModelGrid = true;
    cfg.UseExternalFastKernel = true;
    cfg.ExternalFastKernelByQ = true;
    cfg.ExternalFastKernelFileName = "../../inputs/OutputFiles/sourceKernel_q10_r31_rmax20_coulomb_byQ.root";
    cfg.NumGridBins = 31;
    cfg.gridMinFm = -20.0;
    cfg.gridMaxFm = 20.0;
    cfg.BackgroundNonFemtoQminMeV = 100.0;
    cfg.BackgroundNonFemtoQmaxMeV = 200.0;

    cfg.UseResonanceSource = true;
    cfg.ResonanceWeight1 = 0.6291;
    cfg.ResonanceWeight2 = 0.6291;
    cfg.ResonancePRFileName = "../../inputs/EPOS_data/resonance_PR.root";
    cfg.ResonanceRRFileName = "../../inputs/EPOS_data/resonance_RR.root";
    cfg.ResonanceAlphaInit = 2.0;
    cfg.ResonanceAlphaMin = 1.95;
    cfg.ResonanceAlphaMax = 2.05;
    cfg.FixResonanceAlpha = true;
    cfg.ResonanceUseEllipsoidalCore = true;
    cfg.ResonanceMcIter = 50000;
    cfg.ResonanceRadBins = 51;
    cfg.ResonanceRadMinFm = -40.0;
    cfg.ResonanceRadMaxFm = 40.0;
    cfg.RoutInit = 8.99;
    cfg.RsideInit = 8.00;
    cfg.RlongInit = 9.41;
    if(gaussianControl){
        cfg.UseResonanceSource = false;
        cfg.OutputFileName = "derivedFull_Bminus_cent0to10_anyCPR_resonance6291_q10_r20_weak1557_fit/mt0/qn0/gaussian_control_phi0.root";
    }

    try {
        std::filesystem::create_directories(
            std::filesystem::path(cfg.OutputFileName.Data()).parent_path());
        std::unique_ptr<DLM_CleverReso3D> source = BuildCATS3DResonanceSource(cfg);
        const auto result = FitDataCATS3D(cfg, source.get());
        return result.FitStatus == 0 ? 0 : 3;
    } catch(const std::exception& e) {
        std::cerr << "fit_resonance6291_mt0: " << e.what() << '\n';
        return 1;
    }
}
