#include "FitCATS3D.h"
#include "FitCATS3DRunnerConfig.h"

#include "TROOT.h"

int main()
{
    gROOT->SetBatch(true);
    FitCATS3DConfig cfg = MakeFitCATS3DBaseConfig();
    FitDataCATS3D(cfg);
    return 0;
}
