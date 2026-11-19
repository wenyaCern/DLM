#include "Rtypes.h"
R__LOAD_LIBRARY(/Users/wuwenya/Documents/code/femto_flow/CATS/CATS_install/CATS_3D/build/libFitCATS3D.dylib)

#include "Basics3D.h"

void run_gauss_R20_data_range_check()
{
    Basics_RunGaussianR20DataRangeCheck3D(
        "../OutputFiles/gauss_R20_q1000_check.root",
        1000.0,
        11,
        false,
        true);
}
