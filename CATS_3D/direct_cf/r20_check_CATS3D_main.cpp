#include "Basics3D.h"

#include "TROOT.h"

int main()
{
    gROOT->SetBatch(true);
    Basics_RunGaussianR20DataRangeCheck3D("../../inputs/OutputFiles/gauss_R20_data_range_check.root");
    return 0;
}
