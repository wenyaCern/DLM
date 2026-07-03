#include "Rtypes.h"

R__LOAD_LIBRARY(/Users/wuwenya/Documents/code/femto_flow/CATS/CATS_install/install/lib/libCATS.dylib)
R__LOAD_LIBRARY(libMinuit)
R__LOAD_LIBRARY(/Users/wuwenya/Documents/code/femto_flow/CATS/CATS_install/CATS_3D/build/libFitCATS3D.dylib)

#include "FitCATS3D.h"
#include "FitCATS3DRunnerConfig.h"

void run_fitCATS3D()
{
    FitCATS3DConfig cfg = MakeFitCATS3DBaseConfig();
    FitDataCATS3D(cfg);
}
