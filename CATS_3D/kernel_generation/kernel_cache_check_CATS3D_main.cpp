#include "Basics3D.h"
#include "CATS3DFastKernelCache.h"

#include "CATSconstants.h"

#include "TFile.h"
#include "TH3D.h"

#include <cmath>
#include <iostream>
#include <memory>

namespace {

void ConfigureKitty(CATS3D& kitty,
                    const bool withCoulomb,
                    const double radiusFm)
{
    const unsigned nQ = 10;
    const double qMin = -200.0;
    const double qMax = 200.0;
    kitty.SetMomBinsOut(nQ, qMin, qMax);
    kitty.SetMomBinsSide(nQ, qMin, qMax);
    kitty.SetMomBinsLong(nQ, qMin, qMax);

    const unsigned nR = 11;
    const double rMin = -50.0;
    const double rMax = 50.0;
    kitty.SetGridX(nR, rMin, rMax);
    kitty.SetGridY(nR, rMin, rMax);
    kitty.SetGridZ(nR, rMin, rMax);

    kitty.SetFastKernelMode(true);
    kitty.SetFastKernelSubdivisions(1);

    CATSparameters sourcePars(CATSparameters::tSource3D, 3, true);
    sourcePars.SetParameter(0, radiusFm);
    sourcePars.SetParameter(1, radiusFm);
    sourcePars.SetParameter(2, radiusFm);
    kitty.SetAnaSource(Basics_Source_Gauss3D, sourcePars);

    kitty.SetPdgId(211, 211);
    kitty.SetQuantumStatistics(1);
    kitty.SetQ1Q2(withCoulomb ? 1 : 0);
    kitty.SetRedMass(0.5*Mass_Pic);
    kitty.SetWfMode(withCoulomb ? CATS3D::wfCoulomb : CATS3D::wfPlaneWave);
}

double CompareCats(const CATS3D& ref,
                   const CATS3D& test)
{
    double maxAbs = 0.0;
    for(unsigned ix=0; ix<ref.GetNumMomBinsOut(); ++ix){
        for(unsigned iy=0; iy<ref.GetNumMomBinsSide(); ++iy){
            for(unsigned iz=0; iz<ref.GetNumMomBinsLong(); ++iz){
                const double diff =
                    std::fabs(ref.GetCorrFun(ix, iy, iz) -
                              test.GetCorrFun(ix, iy, iz));
                if(diff > maxAbs) maxAbs = diff;
            }
        }
    }
    return maxAbs;
}

} // namespace

int main(int argc, char** argv)
{
    const bool withCoulomb =
        argc > 1 && TString(argv[1]).Contains("coulomb", TString::kIgnoreCase);
    const double radiusFm = 6.0;

    const TString tag = withCoulomb ? "coulomb" : "noCoulomb";
    const TString outFile = Form("../../inputs/OutputFiles/kernel_cache_check_%s.root", tag.Data());
    const TString kernelFile = Form("../../inputs/OutputFiles/kernel_cache_%s.root", tag.Data());

    CATS3D refKitty;
    ConfigureKitty(refKitty, withCoulomb, radiusFm);
    refKitty.KillTheCat();

    if(!SaveCATS3DFastKernelByQ(refKitty, kernelFile)){
        return 1;
    }

    CATS3D cachedKitty;
    ConfigureKitty(cachedKitty, withCoulomb, radiusFm);
    if(!LoadCATS3DFastKernelByQ(cachedKitty, kernelFile)){
        return 2;
    }
    cachedKitty.KillTheCat();

    const double maxAbs = CompareCats(refKitty, cachedKitty);
    std::cout << "CATS3D kernel cache closure maxAbsDiff=" << maxAbs << "\n";

    std::unique_ptr<TH3D> hRef(Basics_CATS3D_GetCkHisto(
        refKitty, Form("hCk_ref_%s", tag.Data())));
    std::unique_ptr<TH3D> hCached(Basics_CATS3D_GetCkHisto(
        cachedKitty, Form("hCk_cached_%s", tag.Data())));
    std::unique_ptr<TH3D> hDiff(dynamic_cast<TH3D*>(hCached->Clone(
        Form("hCk_diff_%s", tag.Data()))));
    hDiff->Add(hRef.get(), -1.0);

    TFile fout(outFile, "RECREATE");
    hRef->Write();
    hCached->Write();
    hDiff->Write();
    fout.Close();

    std::cout << "Wrote closure file: " << outFile << "\n";
    return maxAbs < 1e-6 ? 0 : 3;
}
