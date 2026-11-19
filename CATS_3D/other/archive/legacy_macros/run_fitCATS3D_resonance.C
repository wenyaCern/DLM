#include "Rtypes.h"

R__LOAD_LIBRARY(/Users/wuwenya/Documents/code/femto_flow/CATS/CATS_install/install/lib/libCATS.dylib)
R__LOAD_LIBRARY(libMinuit)
R__LOAD_LIBRARY(/Users/wuwenya/Documents/code/femto_flow/CATS/CATS_install/CATS_3D/build/libFitCATS3DResonance.dylib)

#include "FitCATS3D.h"
#include "FitCATS3DRunnerConfig.h"
#include "/Users/wuwenya/Documents/code/femto_flow/CATS/CATS_install/install/include/DLM_Source.h"

void run_fitCATS3D_resonance()
{
    FitCATS3DConfig cfg = MakeFitCATS3DResonanceConfig();
    std::unique_ptr<DLM_CleverReso3D> source = BuildCATS3DResonanceSource(cfg);
    FitDataCATS3D(cfg, source.get());
}
