#include "FitCATS3D.h"
#include "FitCATS3DRunnerConfig.h"
#include "Basics3D.h"
#include "/Users/wuwenya/Documents/code/femto_flow/CATS/CATS_install/install/include/CATS3D.h"
#include "/Users/wuwenya/Documents/code/femto_flow/CATS/CATS_install/install/include/CATSconstants.h"
#include "/Users/wuwenya/Documents/code/femto_flow/CATS/CATS_install/install/include/DLM_Source.h"

#include <iostream>
#include <memory>

R__LOAD_LIBRARY(/Users/wuwenya/Documents/code/femto_flow/CATS/CATS_install/CATS_3D/build/libFitCATS3D.dylib)
R__LOAD_LIBRARY(/Users/wuwenya/Documents/code/femto_flow/CATS/CATS_install/CATS_3D/build/libFitCATS3DResonance.dylib)

void profile_resonance_model_test()
{
    auto setupCommon = [](CATS3D& kitty) {
        kitty.SetMomBinsOut(20, -200, 200);
        kitty.SetMomBinsSide(20, -200, 200);
        kitty.SetMomBinsLong(20, -200, 200);
        kitty.SetGridX(15, -20, 20);
        kitty.SetGridY(15, -20, 20);
        kitty.SetGridZ(15, -20, 20);
        kitty.SetFastKernelMode(true);
        kitty.SetFastKernelSubdivisions(2);
        kitty.SetPdgId(211, 211);
        kitty.SetQuantumStatistics(1);
        kitty.SetQ1Q2(0);
        kitty.SetRedMass(0.5 * Mass_Pic);
        kitty.SetWfMode(CATS3D::wfPlaneWave);
    };

    {
        std::cout << "GAUSSIAN SOURCE" << std::endl;
        CATS3D kitty;
        setupCommon(kitty);
        CATSparameters pars(CATSparameters::tSource3D, 3, true);
        pars.SetParameter(0, 5.0);
        pars.SetParameter(1, 5.0);
        pars.SetParameter(2, 5.0);
        kitty.SetAnaSource(Basics_Source_Gauss3D, pars);
        kitty.KillTheCat();
        std::cout << "center=" << kitty.GetCorrFun(10,10,10)
                  << " edge=" << kitty.GetCorrFun(0,0,0) << std::endl;
    }

    {
        std::cout << "RESONANCE SOURCE" << std::endl;
        FitCATS3DConfig cfg = MakeFitCATS3DResonanceConfig();
        cfg.ResonancePRFileName = "./EPOS_data/resonance_PR_PRF.root";
        cfg.ResonanceRRFileName = "./EPOS_data/resonance_RR.root";
        cfg.ResonanceMcIter = 20000;
        cfg.UseResonanceSource = true;
        std::unique_ptr<DLM_CleverReso3D> source = BuildCATS3DResonanceSource(cfg);

        CATS3D kitty;
        setupCommon(kitty);
        kitty.SetAnaSource(CatsSourceForwarder, source.get(), 4);
        kitty.SetAnaSource(0, 5.0);
        kitty.SetAnaSource(1, 5.0);
        kitty.SetAnaSource(2, 5.0);
        kitty.SetAnaSource(3, cfg.ResonanceAlphaInit);
        kitty.KillTheCat();
        std::cout << "center=" << kitty.GetCorrFun(10,10,10)
                  << " edge=" << kitty.GetCorrFun(0,0,0) << std::endl;
    }
}
