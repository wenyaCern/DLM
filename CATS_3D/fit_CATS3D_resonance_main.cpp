#include "FitCATS3D.h"
#include "FitCATS3DRunnerConfig.h"

#include "/Users/wuwenya/Documents/code/femto_flow/CATS/CATS_install/install/include/DLM_Source.h"
#include "TROOT.h"

int main()
{
    gROOT->SetBatch(true);
    FitCATS3DConfig cfg = MakeFitCATS3DResonanceConfig();
    std::unique_ptr<DLM_CleverReso3D> source = BuildCATS3DResonanceSource(cfg);
    FitDataCATS3D(cfg, source.get());
    return 0;
}
