#include "Basics3D.h"
#include "/Users/wuwenya/Documents/code/femto_flow/CATS/CATS_install/install/include/CATS3D.h"
#include "/Users/wuwenya/Documents/code/femto_flow/CATS/CATS_install/install/include/CATSconstants.h"

#include <iostream>

R__LOAD_LIBRARY(/Users/wuwenya/Documents/code/femto_flow/CATS/CATS_install/CATS_3D/build/libFitCATS3D.dylib)

void profile_fast_cache_test()
{
    CATS3D kitty;
    kitty.SetMomBinsOut(20, -200, 200);
    kitty.SetMomBinsSide(20, -200, 200);
    kitty.SetMomBinsLong(20, -200, 200);
    kitty.SetGridX(15, -20, 20);
    kitty.SetGridY(15, -20, 20);
    kitty.SetGridZ(15, -20, 20);
    kitty.SetFastKernelMode(true);
    kitty.SetFastKernelSubdivisions(2);

    CATSparameters pars(CATSparameters::tSource3D, 3, true);
    pars.SetParameter(0, 5.0);
    pars.SetParameter(1, 5.0);
    pars.SetParameter(2, 5.0);
    kitty.SetAnaSource(Basics_Source_Gauss3D, pars);

    kitty.SetPdgId(211, 211);
    kitty.SetQuantumStatistics(1);
    kitty.SetQ1Q2(0);
    kitty.SetRedMass(0.5 * Mass_Pic);
    kitty.SetWfMode(CATS3D::wfPlaneWave);

    std::cout << "FIRST KILL" << std::endl;
    kitty.KillTheCat();

    std::cout << "SECOND KILL" << std::endl;
    kitty.SetAnaSource(0, 5.2, true);
    kitty.SetAnaSource(1, 5.1, true);
    kitty.SetAnaSource(2, 5.3, true);
    kitty.KillTheCat();
}
