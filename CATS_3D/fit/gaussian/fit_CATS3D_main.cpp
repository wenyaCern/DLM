#include "FitCATS3D.h"
#include "FitCATS3DRunnerConfig.h"

#include "TROOT.h"
#include "TFile.h"
#include "TH3D.h"
#include "TNamed.h"
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <memory>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>

namespace {

void BuildBootstrapInput(const FitCATS3DConfig& cfg,
                         const int replica,
                         const int mtIndex = -1,
                         const int nSourceSamples = 21,
                         const int drawsPerReplica = 21,
                         const unsigned baseSeed = 314159u)
{
    std::filesystem::create_directories(
        std::filesystem::path(cfg.InputFileName.Data()).parent_path());
    std::mt19937 rng(baseSeed + static_cast<unsigned>(replica));
    std::uniform_int_distribution<int> pick(0, nSourceSamples-1);
    std::unique_ptr<TH3D> sameSum;
    std::unique_ptr<TH3D> mixedSum;
    std::ostringstream drawList;

    for(int iDraw=0; iDraw<drawsPerReplica; ++iDraw){
        const int iSample = pick(rng);
        if(iDraw) drawList << ',';
        drawList << iSample;
        TFile input(Form("../../inputs/File/subSamp/out_3DCF_mT_cent_qn_phi_tofOnly_samp%d.root",
                         iSample), "READ");
        auto* same = dynamic_cast<TH3D*>(input.Get(cfg.SameHistoName));
        auto* mixed = dynamic_cast<TH3D*>(input.Get(cfg.MixedHistoName));
        if(!same || !mixed){
            throw std::runtime_error("BuildBootstrapInput: missing SAME/MIXED histogram in subSample " +
                                     std::to_string(iSample));
        }
        if(!sameSum){
            sameSum.reset(dynamic_cast<TH3D*>(same->Clone(cfg.SameHistoName)));
            mixedSum.reset(dynamic_cast<TH3D*>(mixed->Clone(cfg.MixedHistoName)));
            sameSum->SetDirectory(nullptr);
            mixedSum->SetDirectory(nullptr);
            sameSum->Reset("ICES");
            mixedSum->Reset("ICES");
        }
        sameSum->Add(same);
        mixedSum->Add(mixed);
    }

    TFile output(cfg.InputFileName, "RECREATE");
    if(output.IsZombie()){
        throw std::runtime_error("BuildBootstrapInput: cannot create " +
                                 std::string(cfg.InputFileName.Data()));
    }
    sameSum->Write(cfg.SameHistoName);
    mixedSum->Write(cfg.MixedHistoName);
    TNamed selection("bootstrap_draw_indices", drawList.str().c_str());
    selection.Write();
    std::cout << "Bootstrap replica " << replica
              << " draws: " << drawList.str() << "\n";
}

} // namespace

int main(int argc, char** argv)
{
    gROOT->SetBatch(true);
    FitCATS3DConfig cfg;
    if(argc >= 2 && std::string(argv[1]) == "bootstrap"){
        if(argc < 3){
            std::cerr << "Usage: ./fit_CATS3D bootstrap REPLICA_NUMBER\n";
            return 2;
        }
        const int replica = std::atoi(argv[2]);
        cfg = MakeFitCATS3DBootstrapConfig(replica);
        BuildBootstrapInput(cfg, replica);
    }
    else if(argc >= 2 && std::string(argv[1]) == "bootstrap-mt"){
        if(argc < 4){
            std::cerr << "Usage: ./fit_CATS3D bootstrap-mt MT_INDEX REPLICA_NUMBER\n";
            return 2;
        }
        const int mtIndex = std::atoi(argv[2]);
        const int replica = std::atoi(argv[3]);
        cfg = MakeFitCATS3DBootstrapMtConfig(mtIndex, replica);
        BuildBootstrapInput(cfg, replica, mtIndex);
    }
    else if(argc >= 2 && std::string(argv[1]) == "bootstrap-mt-nominal"){
        if(argc < 4){
            std::cerr << "Usage: ./fit_CATS3D bootstrap-mt-nominal MT_INDEX REPLICA_NUMBER\n";
            return 2;
        }
        const int mtIndex = std::atoi(argv[2]);
        const int replica = std::atoi(argv[3]);
        // Reuse the already-generated bootstrap SAME/MIXED ROOT file.
        cfg = MakeFitCATS3DNominalBootstrapMtConfig(mtIndex, replica);
    }
    else if(argc >= 2 && std::string(argv[1]) == "full-mt"){
        if(argc < 3){
            std::cerr << "Usage: ./fit_CATS3D full-mt MT_INDEX\n";
            return 2;
        }
        cfg = MakeFitCATS3DMtConfig(std::atoi(argv[2]));
    }
    else if(argc >= 2 && std::string(argv[1]) == "full-mt-oldlambda"){
        if(argc < 3){
            std::cerr << "Usage: ./fit_CATS3D full-mt-oldlambda MT_INDEX\n";
            return 2;
        }
        const int mtIndex = std::atoi(argv[2]);
        cfg = MakeFitCATS3DMtConfig(mtIndex);
        cfg.LambdaFileName =
            "../../inputs/File/lam/lambda_LHC24g3/kStarDepLambda_cent10to30_pion_PbPb.root";
        cfg.OutputFileName = Form(
            "../fit_CATS3D_gaussian_cent1_mt%d_oldLambda_wtCoulomb.root",
            mtIndex);
    }
    else if(argc >= 2 && std::string(argv[1]) == "full-mt-weak14"){
        if(argc < 3){
            std::cerr << "Usage: ./fit_CATS3D full-mt-weak14 MT_INDEX\n";
            return 2;
        }
        const int mtIndex = std::atoi(argv[2]);
        cfg = MakeFitCATS3DMtConfig(mtIndex);
        cfg.OutputFileName = Form(
            "../fit_CATS3D_gaussian_cent1_mt%d_weak14_wtCoulomb.root",
            mtIndex);
    }
    else if(argc >= 2 && std::string(argv[1]) == "full-mt-qmag150to200"){
        if(argc < 3){
            std::cerr << "Usage: ./fit_CATS3D full-mt-qmag150to200 MT_INDEX\n";
            return 2;
        }
        const int mtIndex = std::atoi(argv[2]);
        cfg = MakeFitCATS3DMtConfig(mtIndex);
        cfg.UseSphericalNormalizationRegion = true;
        cfg.normQminMeV = 150.0;
        cfg.normQmaxMeV = 200.0;
        cfg.BackgroundNonFemtoQminMeV = 150.0;
        cfg.OutputFileName = Form(
            "../fit_CATS3D_gaussian_cent1_mt%d_qmag150to200_wtCoulomb.root",
            mtIndex);
    }
    else if(argc >= 2 && std::string(argv[1]) == "full-mt-qmag150to200-radial2"){
        if(argc < 3){
            std::cerr << "Usage: ./fit_CATS3D full-mt-qmag150to200-radial2 MT_INDEX\n";
            return 2;
        }
        const int mtIndex = std::atoi(argv[2]);
        cfg = MakeFitCATS3DMtConfig(mtIndex);
        cfg.UseSphericalNormalizationRegion = true;
        cfg.UseRadialQuadraticBackground = true;
        cfg.normQminMeV = 150.0;
        cfg.normQmaxMeV = 200.0;
        cfg.BackgroundNonFemtoQminMeV = 150.0;
        cfg.OutputFileName = Form(
            "../fit_CATS3D_gaussian_cent1_mt%d_qmag150to200_radial2_wtCoulomb.root",
            mtIndex);
    }
    else if(argc >= 2 && std::string(argv[1]) == "full-mt-qmag120to200-radial2"){
        if(argc < 3){
            std::cerr << "Usage: ./fit_CATS3D full-mt-qmag120to200-radial2 MT_INDEX\n";
            return 2;
        }
        const int mtIndex = std::atoi(argv[2]);
        cfg = MakeFitCATS3DMtConfig(mtIndex);
        cfg.UseSphericalNormalizationRegion = true;
        cfg.UseRadialQuadraticBackground = true;
        cfg.normQminMeV = 120.0;
        cfg.normQmaxMeV = 200.0;
        cfg.BackgroundNonFemtoQminMeV = 120.0;
        cfg.OutputFileName = Form(
            "../fit_CATS3D_gaussian_cent1_mt%d_qmag120to200_radial2_wtCoulomb.root",
            mtIndex);
    }
    else if(argc >= 2 && std::string(argv[1]) == "full-mt-qmag120to200-directional2"){
        if(argc < 3){
            std::cerr << "Usage: ./fit_CATS3D full-mt-qmag120to200-directional2 MT_INDEX\n";
            return 2;
        }
        const int mtIndex = std::atoi(argv[2]);
        cfg = MakeFitCATS3DMtConfig(mtIndex);
        cfg.UseSphericalNormalizationRegion = true;
        cfg.UseDiagonalQuadraticBackground = true;
        cfg.UsePmlBackgroundPrefit = true;
        cfg.normQminMeV = 120.0;
        cfg.normQmaxMeV = 200.0;
        cfg.BackgroundNonFemtoQminMeV = 120.0;
        cfg.OutputFileName = Form(
            "../fit_CATS3D_gaussian_cent1_mt%d_qmag120to200_directional2_wtCoulomb.root",
            mtIndex);
    }
    else if(argc >= 2 && std::string(argv[1]) == "full-mt-qmag120to200-directional2-pml-oldlambda"){
        if(argc < 3){
            std::cerr << "Usage: ./fit_CATS3D full-mt-qmag120to200-directional2-pml-oldlambda MT_INDEX\n";
            return 2;
        }
        const int mtIndex = std::atoi(argv[2]);
        cfg = MakeFitCATS3DMtConfig(mtIndex);
        cfg.UseSphericalNormalizationRegion = true;
        cfg.UseDiagonalQuadraticBackground = true;
        cfg.UsePmlBackgroundPrefit = true;
        cfg.normQminMeV = 120.0;
        cfg.normQmaxMeV = 200.0;
        cfg.BackgroundNonFemtoQminMeV = 120.0;
        cfg.LambdaFileName =
            "../../inputs/File/lam/lambda_LHC24g3/kStarDepLambda_cent10to30_pion_PbPb.root";
        cfg.OutputFileName = Form(
            "../fit_CATS3D_gaussian_cent1_mt%d_qmag120to200_directional2_bkgPML_oldLambda_wtCoulomb.root",
            mtIndex);
    }
    else if(argc >= 2 && std::string(argv[1]) == "full-mt-qmag100to200-directional2"){
        if(argc < 3){
            std::cerr << "Usage: ./fit_CATS3D full-mt-qmag100to200-directional2 MT_INDEX\n";
            return 2;
        }
        const int mtIndex = std::atoi(argv[2]);
        cfg = MakeFitCATS3DMtConfig(mtIndex);
        cfg.UseSphericalNormalizationRegion = true;
        cfg.UseDiagonalQuadraticBackground = true;
        cfg.UsePmlBackgroundPrefit = true;
        cfg.normQminMeV = 100.0;
        cfg.normQmaxMeV = 200.0;
        cfg.BackgroundNonFemtoQminMeV = 100.0;
        cfg.OutputFileName = Form(
            "../fit_CATS3D_gaussian_cent1_mt%d_qmag100to200_directional2_wtCoulomb.root",
            mtIndex);
    }
    else if(argc >= 2 && std::string(argv[1]) == "full-mt-qmag150to200-directional2-pml"){
        if(argc < 3){
            std::cerr << "Usage: ./fit_CATS3D full-mt-qmag150to200-directional2-pml MT_INDEX\n";
            return 2;
        }
        const int mtIndex = std::atoi(argv[2]);
        cfg = MakeFitCATS3DMtConfig(mtIndex);
        cfg.UseSphericalNormalizationRegion = true;
        cfg.UseDiagonalQuadraticBackground = true;
        cfg.UsePmlBackgroundPrefit = true;
        cfg.normQminMeV = 150.0;
        cfg.normQmaxMeV = 200.0;
        cfg.BackgroundNonFemtoQminMeV = 150.0;
        cfg.OutputFileName = Form(
            "../fit_CATS3D_gaussian_cent1_mt%d_qmag150to200_directional2_bkgPML_wtCoulomb.root",
            mtIndex);
    }
    else if(argc >= 2 && std::string(argv[1]) == "full-mt-bkgqmag120to200-directional4-pml-oldlambda"){
        if(argc < 3){
            std::cerr << "Usage: ./fit_CATS3D full-mt-bkgqmag120to200-directional4-pml-oldlambda MT_INDEX\n";
            return 2;
        }
        const int mtIndex = std::atoi(argv[2]);
        cfg = MakeFitCATS3DMtConfig(mtIndex);
        // Keep the established normalization configuration untouched.
        cfg.UseDiagonalQuarticBackground = true;
        cfg.UseSphericalBackgroundRegion = true;
        cfg.UsePmlBackgroundPrefit = true;
        cfg.BackgroundNonFemtoQminMeV = 120.0;
        cfg.BackgroundNonFemtoQmaxMeV = 200.0;
        cfg.LambdaFileName =
            "../../inputs/File/lam/lambda_LHC24g3/kStarDepLambda_cent10to30_pion_PbPb.root";
        cfg.OutputFileName = Form(
            "../fit_CATS3D_gaussian_cent1_mt%d_bkgqmag120to200_directional4_bkgPML_oldLambda_wtCoulomb.root",
            mtIndex);
    }
    else if(argc >= 2 && std::string(argv[1]) == "full-mt-bkghybrid100to200-directional2-pml-oldlambda"){
        if(argc < 3){
            std::cerr << "Usage: ./fit_CATS3D full-mt-bkghybrid100to200-directional2-pml-oldlambda MT_INDEX\n";
            return 2;
        }
        const int mtIndex = std::atoi(argv[2]);
        cfg = MakeFitCATS3DMtConfig(mtIndex);
        // Normalization stays exactly as configured by MakeFitCATS3DMtConfig.
        cfg.UseDiagonalQuadraticBackground = true;
        cfg.UseHybridBackgroundRegion = true;
        cfg.UsePmlBackgroundPrefit = true;
        cfg.BackgroundNonFemtoQminMeV = 100.0;
        cfg.BackgroundNonFemtoQmaxMeV = 200.0;
        cfg.sliceCutAbs = 0.04;
        cfg.LambdaFileName =
            "../../inputs/File/lam/lambda_LHC24g3/kStarDepLambda_cent10to30_pion_PbPb.root";
        cfg.weakDecayFraction = 0.12;
        cfg.OutputFileName = Form(
            "../fit_CATS3D_gaussian_cent1_mt%d_bkghybrid100to200_directional2_bkgPML_oldLambda_wtCoulomb.root",
            mtIndex);
    }
    else if(argc >= 2 && std::string(argv[1]) == "sub"){
        if(argc < 3){
            std::cerr << "Usage: ./fit_CATS3D sub SAMPLE_NUMBER\n";
            return 2;
        }
        cfg = MakeFitCATS3DSubSampleConfig(std::atoi(argv[2]));
    }
    else {
        cfg = MakeFitCATS3DBaseConfig();
    }
    FitDataCATS3D(cfg);
    return 0;
}
