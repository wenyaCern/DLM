#include "Basics3D.h"
#include "CATS3DFastKernelCache.h"
#include "FitCATS3D.h"
#include "FitCATS3DRunnerConfig.h"

#include "/Users/wuwenya/Documents/code/femto_flow/CATS/CATS_install/install/include/CATSconstants.h"

#include "TH3D.h"
#include "TString.h"

#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <memory>

namespace {

void ConfigureKernelKitty(CATS3D& kitty,
                          const TH3D& data,
                          const FitCATS3DConfig& cfg)
{
    kitty.SetMomBinsOut(data.GetNbinsX(),
                        data.GetXaxis()->GetXmin()*cfg.qAxisScaleToMeV,
                        data.GetXaxis()->GetXmax()*cfg.qAxisScaleToMeV);
    kitty.SetMomBinsSide(data.GetNbinsY(),
                         data.GetYaxis()->GetXmin()*cfg.qAxisScaleToMeV,
                         data.GetYaxis()->GetXmax()*cfg.qAxisScaleToMeV);
    kitty.SetMomBinsLong(data.GetNbinsZ(),
                         data.GetZaxis()->GetXmin()*cfg.qAxisScaleToMeV,
                         data.GetZaxis()->GetXmax()*cfg.qAxisScaleToMeV);

    kitty.SetGridX(cfg.NumGridBins, cfg.gridMinFm, cfg.gridMaxFm);
    kitty.SetGridY(cfg.NumGridBins, cfg.gridMinFm, cfg.gridMaxFm);
    kitty.SetGridZ(cfg.NumGridBins, cfg.gridMinFm, cfg.gridMaxFm);

    kitty.SetFastKernelMode(true);
    kitty.SetFastKernelSubdivisions(cfg.FastKernelSubdivisions);

    CATSparameters sourcePars(CATSparameters::tSource3D, 3, true);
    sourcePars.SetParameter(0, cfg.RoutInit);
    sourcePars.SetParameter(1, cfg.RsideInit);
    sourcePars.SetParameter(2, cfg.RlongInit);
    kitty.SetAnaSource(Basics_Source_Gauss3D, sourcePars);

    if(cfg.Identical) kitty.SetPdgId(211, 211);
    else              kitty.SetPdgId(-211, 211);

    kitty.SetQuantumStatistics(cfg.Identical ? 1 : 0);
    kitty.SetQ1Q2(cfg.WithCoulomb ? 1 : 0);
    kitty.SetRedMass(0.5*Mass_Pic);
    kitty.SetWfMode(cfg.WithCoulomb ? CATS3D::wfCoulomb : CATS3D::wfPlaneWave);
}

} // namespace

int main(int argc, char** argv)
{
    FitCATS3DConfig cfg = MakeFitCATS3DBaseConfig();

    cfg.RebinFactor = 1;
    cfg.UseFastKernel = true;
    cfg.FastKernelSubdivisions = 1;
    cfg.NumGridBins = 31;
    cfg.gridMinFm = -50.0;
    cfg.gridMaxFm = 50.0;
    cfg.WithCoulomb = false;

    if(argc > 1){
        const TString mode(argv[1]);
        cfg.WithCoulomb =
            mode.EqualTo("coulomb", TString::kIgnoreCase) ||
            mode.EqualTo("withCoulomb", TString::kIgnoreCase) ||
            mode.EqualTo("on", TString::kIgnoreCase);
    }
    if(argc > 2){
        cfg.NumGridBins = unsigned(std::max(1, std::atoi(argv[2])));
    }

    const TString tag = cfg.WithCoulomb ? "coulomb" : "noCoulomb";
    const TString outFile = Form("../OutputFiles/dataKernel_qData_r%u_%s_byQ.root",
                                 cfg.NumGridBins, tag.Data());

    std::unique_ptr<TH3D> data(BuildDataCk3D(cfg));
    if(!data){
        std::cout << "build_data_kernel_CATS3D: failed to build data template.\n";
        return 1;
    }

    CATS3D kitty;
    ConfigureKernelKitty(kitty, *data, cfg);

    const size_t nQ = kitty.GetFastKernelMomSize();
    const size_t nR = kitty.GetFastKernelGridSize();
    const double approxGB = double(nQ*nR*sizeof(float))/1024.0/1024.0/1024.0;
    std::cout << "Building CATS3D by-q kernel"
              << " qBins=(" << kitty.GetNumMomBinsOut()
              << "," << kitty.GetNumMomBinsSide()
              << "," << kitty.GetNumMomBinsLong() << ")"
              << " rBins=(" << kitty.GetNumGridBinsX()
              << "," << kitty.GetNumGridBinsY()
              << "," << kitty.GetNumGridBinsZ() << ")"
              << " rRange=[" << cfg.gridMinFm << "," << cfg.gridMaxFm << "] fm"
              << " subdivisions=" << cfg.FastKernelSubdivisions
              << " withCoulomb=" << cfg.WithCoulomb
              << " totalFloats=" << nQ*nR
              << " rawSize~" << approxGB << " GiB\n";

    if(!SaveCATS3DFastKernelByQ(kitty, outFile)){
        return 2;
    }

    CATS3D checkKitty;
    ConfigureKernelKitty(checkKitty, *data, cfg);
    if(!LoadCATS3DFastKernelByQ(checkKitty, outFile)){
        return 3;
    }

    std::cout << "Kernel workflow check loaded successfully: " << outFile << "\n";
    return 0;
}
