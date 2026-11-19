#include "Basics3D.h"
#include "CATS3DFastKernelCache.h"

#include "CATSconstants.h"

#include <iostream>

int main() {
    constexpr unsigned kMomBins = 20;
    constexpr unsigned kSourceBins = 31;
    constexpr double kQMinMeV = -100.0;
    constexpr double kQMaxMeV = 100.0;
    constexpr double kRMinFm = -20.0;
    constexpr double kRMaxFm = 20.0;
    const TString output = "../../inputs/OutputFiles/sourceKernel_q10_r31_rmax20_coulomb_byQ.root";

    CATS3D kitty;
    kitty.SetMomBinsOut(kMomBins, kQMinMeV, kQMaxMeV);
    kitty.SetMomBinsSide(kMomBins, kQMinMeV, kQMaxMeV);
    kitty.SetMomBinsLong(kMomBins, kQMinMeV, kQMaxMeV);
    kitty.SetGridX(kSourceBins, kRMinFm, kRMaxFm);
    kitty.SetGridY(kSourceBins, kRMinFm, kRMaxFm);
    kitty.SetGridZ(kSourceBins, kRMinFm, kRMaxFm);
    kitty.SetFastKernelMode(true);
    kitty.SetFastKernelSubdivisions(1);

    // The wave-function kernel is independent of these source radii; they
    // provide a valid source configuration for the CATS3D build.
    CATSparameters source(CATSparameters::tSource3D, 3, true);
    source.SetParameter(0, 8.0);
    source.SetParameter(1, 8.0);
    source.SetParameter(2, 8.0);
    kitty.SetAnaSource(Basics_Source_Gauss3D, source);
    kitty.SetPdgId(211, 211);
    kitty.SetQuantumStatistics(1);
    kitty.SetQ1Q2(1);
    kitty.SetRedMass(0.5 * Mass_Pic);
    kitty.SetWfMode(CATS3D::wfCoulomb);

    std::cout << "Building " << output << " with qBins=" << kMomBins
              << "^3 qRange=[" << kQMinMeV << ',' << kQMaxMeV
              << "] MeV/c rBins=" << kSourceBins << "^3 rRange=["
              << kRMinFm << ',' << kRMaxFm << "] fm\n" << std::flush;
    return SaveCATS3DFastKernelByQ(kitty, output) ? 0 : 1;
}
