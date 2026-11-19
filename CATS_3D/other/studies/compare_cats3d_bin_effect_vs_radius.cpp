#include "Basics3D.h"
#include "/Users/wuwenya/Documents/code/femto_flow/CATS/CATS_install/install/include/CATSconstants.h"

#include <iomanip>
#include <iostream>

namespace {
constexpr double kLambda = 0.83 * (1.0 - 0.14) * (1.0 - 0.14);

void Configure(CATS3D& kitty, const double radius, const unsigned qBins) {
    kitty.SetMomBinsOut(qBins, 0.0, 40.0);
    kitty.SetMomBinsSide(qBins, 0.0, 40.0);
    kitty.SetMomBinsLong(qBins, 0.0, 40.0);
    kitty.SetGridX(41, -40.0, 40.0);
    kitty.SetGridY(41, -40.0, 40.0);
    kitty.SetGridZ(41, -40.0, 40.0);
    kitty.SetFastKernelMode(true);
    kitty.SetFastKernelSubdivisions(1);
    CATSparameters source(CATSparameters::tSource3D, 3, true);
    for (unsigned axis = 0; axis < 3; ++axis) source.SetParameter(axis, radius);
    kitty.SetAnaSource(Basics_Source_Gauss3D, source);
    kitty.SetPdgId(211, 211);
    kitty.SetQuantumStatistics(1);
    kitty.SetQ1Q2(1);
    kitty.SetRedMass(0.5 * Mass_Pic);
    kitty.SetWfMode(CATS3D::wfCoulomb);
}

double C(const CATS3D& kitty, unsigned ix, unsigned iy, unsigned iz) {
    return 1.0 + kLambda * (kitty.GetCorrFun(ix, iy, iz) - 1.0);
}

void Compare(const double radius) {
    CATS3D coarse;
    Configure(coarse, radius, 2); // 20 MeV/c bins in [0,40]^3
    coarse.KillTheCat();
    CATS3D fine;
    Configure(fine, radius, 4); // 10 MeV/c bins in [0,40]^3
    fine.KillTheCat();

    const double center20 = C(coarse, 0, 0, 0);
    double average10 = 0.0;
    for (unsigned ix = 0; ix < 2; ++ix)
      for (unsigned iy = 0; iy < 2; ++iy)
        for (unsigned iz = 0; iz < 2; ++iz)
          average10 += C(fine, ix, iy, iz) / 8.0;

    std::cout << radius << ',' << center20 << ',' << average10 << ','
              << average10 - center20 << ','
              << (average10 - center20) / (center20 - 1.0) << '\n';
}
}

int main() {
    std::cout << std::setprecision(9);
    std::cout << "R_fm,C20_center,C10_average,delta,delta_over_C20_minus_1\n";
    Compare(8.0);
    Compare(3.0);
}
