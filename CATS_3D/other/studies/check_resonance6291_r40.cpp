#include "Basics3D.h"
#include "FitCATS3D.h"
#include "FitCATS3DRunnerConfig.h"

#include "/Users/wuwenya/Documents/code/femto_flow/CATS/CATS_install/install/include/CATSconstants.h"
#include "/Users/wuwenya/Documents/code/femto_flow/CATS/CATS_install/install/include/DLM_Source.h"

#include <fstream>
#include <iostream>
#include <memory>
#include <string>

int main(int argc, char** argv)
{
    const bool gaussian = argc == 2 && std::string(argv[1]) == "gaussian";
    auto cfg = MakeFitCATS3DDerivedFullBminusConfig(0, 0, 0);
    cfg.UseResonanceSource = true;
    cfg.ResonanceWeight1 = cfg.ResonanceWeight2 = 0.6291;
    cfg.ResonancePRFileName = "EPOS_data/resonance_PR.root";
    cfg.ResonanceRRFileName = "EPOS_data/resonance_RR.root";
    cfg.ResonanceMcIter = 50000;
    cfg.ResonanceRadBins = 51;
    cfg.ResonanceRadMinFm = -40.0;
    cfg.ResonanceRadMaxFm = 40.0;
    cfg.ResonanceAlphaInit = 2.0;
    cfg.ResonanceAlphaMin = 1.95;
    cfg.ResonanceAlphaMax = 2.05;
    cfg.RoutInit = 9.1056;
    cfg.RsideInit = 7.26473;
    cfg.RlongInit = 8.95836;
    if(gaussian) cfg.UseResonanceSource = false;
    std::unique_ptr<DLM_CleverReso3D> source = BuildCATS3DResonanceSource(cfg);
    const double lambda = 0.83 * (1.0 - 0.1557) * (1.0 - 0.1557);
    const char* label[] = {"out", "side", "long"};
    const TString output = gaussian
        ? "derivedFull_Bminus_cent0to10_anyCPR_resonance6291_q10_r20_weak1557_fit/mt0/qn0/gaussian_r20_vs_r40_selected_q.csv"
        : "derivedFull_Bminus_cent0to10_anyCPR_resonance6291_q10_r20_weak1557_fit/mt0/qn0/r20_vs_r40_selected_q.csv";
    std::ofstream csv(output.Data());
    csv << "rmax_fm,axis,qcenter_GeV,modelC\n";

    for(const double rmax : {20.0, 40.0}) {
        CATS3D kitty;
        kitty.SetMomBinsOut(4, 0.0, 40.0);
        kitty.SetMomBinsSide(4, 0.0, 40.0);
        kitty.SetMomBinsLong(4, 0.0, 40.0);
        kitty.SetGridX(31, -rmax, rmax);
        kitty.SetGridY(31, -rmax, rmax);
        kitty.SetGridZ(31, -rmax, rmax);
        kitty.SetFastKernelMode(true);
        kitty.SetFastKernelSubdivisions(1);
        if(gaussian){
            CATSparameters pars(CATSparameters::tSource3D, 3, true);
            pars.SetParameter(0, cfg.RoutInit);
            pars.SetParameter(1, cfg.RsideInit);
            pars.SetParameter(2, cfg.RlongInit);
            kitty.SetAnaSource(Basics_Source_Gauss3D, pars);
        } else {
            kitty.SetAnaSource(CatsSourceForwarder, source.get(), 4);
            kitty.SetAnaSource(0, cfg.RoutInit);
            kitty.SetAnaSource(1, cfg.RsideInit);
            kitty.SetAnaSource(2, cfg.RlongInit);
            kitty.SetAnaSource(3, 2.0);
        }
        kitty.SetPdgId(211, 211);
        kitty.SetQuantumStatistics(1);
        kitty.SetQ1Q2(1);
        kitty.SetRedMass(0.5 * Mass_Pic);
        kitty.SetWfMode(CATS3D::wfCoulomb);
        kitty.KillTheCat();
        for(unsigned axis=0; axis<3; ++axis){
            for(unsigned shoulder=0; shoulder<2; ++shoulder){
                double mean=0.0;
                for(unsigned i=0;i<2;++i)
                    for(unsigned j=0;j<2;++j)
                        for(unsigned k=0;k<2;++k){
                            unsigned idx[3]={i,j,k};
                            idx[axis]+=2*shoulder;
                            mean+=kitty.GetCorrFun(idx[0],idx[1],idx[2])/8.0;
                        }
                const double model=1.0+lambda*(mean-1.0);
                csv<<rmax<<','<<label[axis]<<','<<(shoulder?0.03:0.01)<<','<<model<<'\n';
                std::cout<<rmax<<' '<<label[axis]<<' '<<(shoulder?0.03:0.01)
                         <<' '<<model<<'\n';
            }
        }
    }
}
