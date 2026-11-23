#include "FitCATS3DRunnerConfig.h"

#include "TROOT.h"

#include <filesystem>
#include <array>
#include <iostream>
#include <stdexcept>
#include <vector>

int main(int argc, char** argv)
{
    gROOT->SetBatch(true);
    if(argc < 4){
        std::cerr << "Usage: fit_CATS3D_derived_full MT QN PHI [options], or --batch MT QN_FIRST QN_LAST PHI_FIRST PHI_LAST [options]\n";
        return 2;
    }
    try{
        const bool batch=std::string(argv[1])=="--batch";
        std::vector<std::array<int,3>> cells;
        int firstOption=4;
        if(batch){
            if(argc<7) throw std::runtime_error("Batch needs MT QN_FIRST QN_LAST PHI_FIRST PHI_LAST");
            const int mt=std::stoi(argv[2]);
            const int qnFirst=std::stoi(argv[3]), qnLast=std::stoi(argv[4]);
            const int phiFirst=std::stoi(argv[5]), phiLast=std::stoi(argv[6]);
            if(mt<0 || mt>2 || qnFirst<0 || qnLast>4 || qnFirst>qnLast ||
               phiFirst<0 || phiLast>5 || phiFirst>phiLast)
                throw std::runtime_error("Invalid batch cell range");
            for(int qn=qnFirst;qn<=qnLast;++qn)
                for(int phi=phiFirst;phi<=phiLast;++phi)
                    cells.push_back({mt,qn,phi});
            firstOption=7;
        } else {
            cells.push_back({std::stoi(argv[1]),std::stoi(argv[2]),std::stoi(argv[3])});
        }
        int failures=0;
        for(const auto& cell:cells){
        auto cfg = MakeFitCATS3DDerivedFullBminusConfig(cell[0],cell[1],cell[2]);
        const TString batchOutput=cfg.OutputFileName;
        for(int i=firstOption; i<argc; ++i){
            const std::string arg(argv[i]);
            if(arg.rfind("adaptive-kernel=",0) == 0){
                cfg.UseAdaptiveKernel=true;
                cfg.AdaptiveKernelFileName=arg.substr(16).c_str();
                cfg.UseFastKernel=false;
                cfg.UseExternalFastKernel=false;
                cfg.UseFineModelGrid=false;
                cfg.UseModelCache=true;
                cfg.ModelCacheToleranceFm=0.02;
            }
            else if(arg == "lambda-after-dca-cut"){
                cfg.LambdaHistoName = "lambdaAfterDCACut";
            }
            else if(arg.rfind("weakfraction=",0) == 0){
                const double fraction = std::stod(arg.substr(13));
                if(fraction < 0.0 || fraction >= 1.0)
                    throw std::runtime_error("weakfraction must be within [0, 1)");
                cfg.weakDecayFraction = fraction;
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
            else if(arg == "no-simplex"){
                cfg.RunMinuitSimplexFirst = false;
            }
            else if(arg == "kernel-fitbox-only"){
                cfg.OutputBaselineOutsideFitBox = true;
            }
            else if(arg.rfind("output=",0) == 0){
                if(batch) throw std::runtime_error("Use CATS3D_FIT_OUTPUT_DIR for batch output");
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
        if(batch) cfg.OutputFileName=batchOutput;
        std::filesystem::create_directories(
            std::filesystem::path(cfg.OutputFileName.Data()).parent_path());
        std::cout << "\n=== CATS3D cell mT=" << cell[0] << " qn=" << cell[1]
                  << " phi=" << cell[2] << " ===\n";
        try{
            const auto result = FitDataCATS3D(cfg);
            if(result.FitStatus!=0) ++failures;
        }catch(const std::exception& e){
            std::cerr << "Cell failed: " << e.what() << '\n';
            ++failures;
            if(!batch) throw;
        }
        }
        std::cout << "CATS3D completed " << cells.size()-failures << "/"
                  << cells.size() << " cells\n";
        return failures ? 3 : 0;
    }
    catch(const std::exception& e){
        std::cerr << "fit_CATS3D_derived_full: " << e.what() << '\n';
        return 1;
    }
}
