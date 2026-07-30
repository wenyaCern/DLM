#include "CATS3DFastKernelCache.h"

#include "TFile.h"
#include "TTree.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <vector>

namespace {

constexpr double kAbsTol = 1e-9;
constexpr double kRelTol = 1e-9;

bool NearlyEqual(const double a, const double b)
{
    const double scale = std::max({1.0, std::fabs(a), std::fabs(b)});
    return std::fabs(a-b) <= kAbsTol + kRelTol*scale;
}

std::vector<double> MomentumEdgesOut(const CATS3D& kitty)
{
    std::vector<double> edges(kitty.GetNumMomBinsOut()+1, 0.0);
    for(unsigned i=0; i<kitty.GetNumMomBinsOut(); ++i){
        edges[i] = kitty.GetMomBinLowEdgeOut(i);
    }
    if(kitty.GetNumMomBinsOut()){
        edges.back() = kitty.GetMomBinUpEdgeOut(kitty.GetNumMomBinsOut()-1);
    }
    return edges;
}

std::vector<double> MomentumEdgesSide(const CATS3D& kitty)
{
    std::vector<double> edges(kitty.GetNumMomBinsSide()+1, 0.0);
    for(unsigned i=0; i<kitty.GetNumMomBinsSide(); ++i){
        edges[i] = kitty.GetMomBinLowEdgeSide(i);
    }
    if(kitty.GetNumMomBinsSide()){
        edges.back() = kitty.GetMomBinUpEdgeSide(kitty.GetNumMomBinsSide()-1);
    }
    return edges;
}

std::vector<double> MomentumEdgesLong(const CATS3D& kitty)
{
    std::vector<double> edges(kitty.GetNumMomBinsLong()+1, 0.0);
    for(unsigned i=0; i<kitty.GetNumMomBinsLong(); ++i){
        edges[i] = kitty.GetMomBinLowEdgeLong(i);
    }
    if(kitty.GetNumMomBinsLong()){
        edges.back() = kitty.GetMomBinUpEdgeLong(kitty.GetNumMomBinsLong()-1);
    }
    return edges;
}

std::vector<double> GridCentersX(const CATS3D& kitty)
{
    std::vector<double> centers(kitty.GetNumGridBinsX(), 0.0);
    for(unsigned i=0; i<kitty.GetNumGridBinsX(); ++i){
        centers[i] = kitty.GetGridCenterX(i);
    }
    return centers;
}

std::vector<double> GridCentersY(const CATS3D& kitty)
{
    std::vector<double> centers(kitty.GetNumGridBinsY(), 0.0);
    for(unsigned i=0; i<kitty.GetNumGridBinsY(); ++i){
        centers[i] = kitty.GetGridCenterY(i);
    }
    return centers;
}

std::vector<double> GridCentersZ(const CATS3D& kitty)
{
    std::vector<double> centers(kitty.GetNumGridBinsZ(), 0.0);
    for(unsigned i=0; i<kitty.GetNumGridBinsZ(); ++i){
        centers[i] = kitty.GetGridCenterZ(i);
    }
    return centers;
}

bool SameVector(const std::vector<double>& a,
                const std::vector<double>& b,
                const char* label)
{
    if(a.size() != b.size()){
        std::cout << "CATS3D fast-kernel cache mismatch: " << label
                  << " size " << a.size() << " vs " << b.size() << "\n";
        return false;
    }
    for(size_t i=0; i<a.size(); ++i){
        if(!NearlyEqual(a[i], b[i])){
            std::cout << "CATS3D fast-kernel cache mismatch: " << label
                      << "[" << i << "] " << a[i] << " vs " << b[i] << "\n";
            return false;
        }
    }
    return true;
}

bool CheckScalar(const char* label,
                 const double stored,
                 const double current)
{
    if(NearlyEqual(stored, current)) return true;
    std::cout << "CATS3D fast-kernel cache mismatch: " << label
              << " " << stored << " vs " << current << "\n";
    return false;
}

bool CheckScalar(const char* label,
                 const unsigned stored,
                 const unsigned current)
{
    if(stored == current) return true;
    std::cout << "CATS3D fast-kernel cache mismatch: " << label
              << " " << stored << " vs " << current << "\n";
    return false;
}

} // namespace

bool SaveCATS3DFastKernel(CATS3D& kitty,
                          const TString& fileName,
                          const TString& treeName,
                          const bool buildIfNeeded)
{
    if(buildIfNeeded && !kitty.IsFastKernelReady()){
        kitty.BuildFastKernelOnly();
    }
    if(!kitty.IsFastKernelReady()){
        std::cout << "SaveCATS3DFastKernel: fast kernel is not ready.\n";
        return false;
    }

    const std::vector<float>& kernel = kitty.GetFastWaveFunctionKernel();
    if(kernel.size() != kitty.GetFastKernelStorageSize()){
        std::cout << "SaveCATS3DFastKernel: kernel size mismatch before write.\n";
        return false;
    }

    std::unique_ptr<TFile> fout(TFile::Open(fileName, "RECREATE"));
    if(!fout || fout->IsZombie()){
        std::cout << "SaveCATS3DFastKernel: cannot open " << fileName << "\n";
        return false;
    }

    unsigned version = 1;
    unsigned subdivisions = kitty.GetFastKernelSubdivisions();
    unsigned nQOut = kitty.GetNumMomBinsOut();
    unsigned nQSide = kitty.GetNumMomBinsSide();
    unsigned nQLong = kitty.GetNumMomBinsLong();
    unsigned nRX = kitty.GetNumGridBinsX();
    unsigned nRY = kitty.GetNumGridBinsY();
    unsigned nRZ = kitty.GetNumGridBinsZ();
    unsigned maxPw = kitty.GetMaxPw();
    unsigned wfMode = unsigned(kitty.GetWfMode());
    unsigned identical = kitty.GetIdenticalParticles() ? 1u : 0u;
    int q1q2 = kitty.GetQ1Q2();
    double redMass = kitty.GetRedMass();
    double gridStepX = kitty.GetGridStepX();
    double gridStepY = kitty.GetGridStepY();
    double gridStepZ = kitty.GetGridStepZ();

    std::vector<double> qOutEdges = MomentumEdgesOut(kitty);
    std::vector<double> qSideEdges = MomentumEdgesSide(kitty);
    std::vector<double> qLongEdges = MomentumEdgesLong(kitty);
    std::vector<double> rOutCenters = GridCentersX(kitty);
    std::vector<double> rSideCenters = GridCentersY(kitty);
    std::vector<double> rLongCenters = GridCentersZ(kitty);
    std::vector<float> kernelCopy = kernel;

    TTree tree(treeName, treeName);
    tree.Branch("version", &version);
    tree.Branch("subdivisions", &subdivisions);
    tree.Branch("nQOut", &nQOut);
    tree.Branch("nQSide", &nQSide);
    tree.Branch("nQLong", &nQLong);
    tree.Branch("nRX", &nRX);
    tree.Branch("nRY", &nRY);
    tree.Branch("nRZ", &nRZ);
    tree.Branch("maxPw", &maxPw);
    tree.Branch("wfMode", &wfMode);
    tree.Branch("identical", &identical);
    tree.Branch("q1q2", &q1q2);
    tree.Branch("redMass", &redMass);
    tree.Branch("gridStepX", &gridStepX);
    tree.Branch("gridStepY", &gridStepY);
    tree.Branch("gridStepZ", &gridStepZ);
    tree.Branch("qOutEdges", &qOutEdges);
    tree.Branch("qSideEdges", &qSideEdges);
    tree.Branch("qLongEdges", &qLongEdges);
    tree.Branch("rOutCenters", &rOutCenters);
    tree.Branch("rSideCenters", &rSideCenters);
    tree.Branch("rLongCenters", &rLongCenters);
    tree.Branch("kernel", &kernelCopy);
    tree.Fill();
    tree.Write("", TObject::kOverwrite);
    fout->Close();

    std::cout << "Saved CATS3D fast kernel: " << fileName
              << " entries=" << kernelCopy.size()
              << " qBins=(" << nQOut << "," << nQSide << "," << nQLong << ")"
              << " rBins=(" << nRX << "," << nRY << "," << nRZ << ")"
              << " subdivisions=" << subdivisions << "\n";
    return true;
}

bool LoadCATS3DFastKernel(CATS3D& kitty,
                          const TString& fileName,
                          const TString& treeName)
{
    std::unique_ptr<TFile> fin(TFile::Open(fileName, "READ"));
    if(!fin || fin->IsZombie()){
        std::cout << "LoadCATS3DFastKernel: cannot open " << fileName << "\n";
        return false;
    }

    TTree* tree = dynamic_cast<TTree*>(fin->Get(treeName));
    if(!tree){
        std::cout << "LoadCATS3DFastKernel: cannot find tree " << treeName
                  << " in " << fileName << "\n";
        return false;
    }
    if(tree->GetEntries() <= 0){
        std::cout << "LoadCATS3DFastKernel: tree is empty.\n";
        return false;
    }

    unsigned version = 0;
    unsigned subdivisions = 0;
    unsigned nQOut = 0;
    unsigned nQSide = 0;
    unsigned nQLong = 0;
    unsigned nRX = 0;
    unsigned nRY = 0;
    unsigned nRZ = 0;
    unsigned maxPw = 0;
    unsigned wfMode = 0;
    unsigned identical = 0;
    int q1q2 = 0;
    double redMass = 0.0;
    double gridStepX = 0.0;
    double gridStepY = 0.0;
    double gridStepZ = 0.0;
    std::vector<double>* qOutEdges = nullptr;
    std::vector<double>* qSideEdges = nullptr;
    std::vector<double>* qLongEdges = nullptr;
    std::vector<double>* rOutCenters = nullptr;
    std::vector<double>* rSideCenters = nullptr;
    std::vector<double>* rLongCenters = nullptr;
    std::vector<float>* kernel = nullptr;

    tree->SetBranchAddress("version", &version);
    tree->SetBranchAddress("subdivisions", &subdivisions);
    tree->SetBranchAddress("nQOut", &nQOut);
    tree->SetBranchAddress("nQSide", &nQSide);
    tree->SetBranchAddress("nQLong", &nQLong);
    tree->SetBranchAddress("nRX", &nRX);
    tree->SetBranchAddress("nRY", &nRY);
    tree->SetBranchAddress("nRZ", &nRZ);
    tree->SetBranchAddress("maxPw", &maxPw);
    tree->SetBranchAddress("wfMode", &wfMode);
    tree->SetBranchAddress("identical", &identical);
    tree->SetBranchAddress("q1q2", &q1q2);
    tree->SetBranchAddress("redMass", &redMass);
    tree->SetBranchAddress("gridStepX", &gridStepX);
    tree->SetBranchAddress("gridStepY", &gridStepY);
    tree->SetBranchAddress("gridStepZ", &gridStepZ);
    tree->SetBranchAddress("qOutEdges", &qOutEdges);
    tree->SetBranchAddress("qSideEdges", &qSideEdges);
    tree->SetBranchAddress("qLongEdges", &qLongEdges);
    tree->SetBranchAddress("rOutCenters", &rOutCenters);
    tree->SetBranchAddress("rSideCenters", &rSideCenters);
    tree->SetBranchAddress("rLongCenters", &rLongCenters);
    tree->SetBranchAddress("kernel", &kernel);
    tree->GetEntry(0);

    bool ok = true;
    ok = CheckScalar("version", version, 1u) && ok;
    ok = CheckScalar("subdivisions", subdivisions, kitty.GetFastKernelSubdivisions()) && ok;
    ok = CheckScalar("nQOut", nQOut, kitty.GetNumMomBinsOut()) && ok;
    ok = CheckScalar("nQSide", nQSide, kitty.GetNumMomBinsSide()) && ok;
    ok = CheckScalar("nQLong", nQLong, kitty.GetNumMomBinsLong()) && ok;
    ok = CheckScalar("nRX", nRX, kitty.GetNumGridBinsX()) && ok;
    ok = CheckScalar("nRY", nRY, kitty.GetNumGridBinsY()) && ok;
    ok = CheckScalar("nRZ", nRZ, kitty.GetNumGridBinsZ()) && ok;
    ok = CheckScalar("maxPw", maxPw, unsigned(kitty.GetMaxPw())) && ok;
    ok = CheckScalar("wfMode", wfMode, unsigned(kitty.GetWfMode())) && ok;
    ok = CheckScalar("identical", identical, kitty.GetIdenticalParticles() ? 1u : 0u) && ok;
    if(q1q2 != kitty.GetQ1Q2()){
        std::cout << "CATS3D fast-kernel cache mismatch: q1q2 "
                  << q1q2 << " vs " << kitty.GetQ1Q2() << "\n";
        ok = false;
    }
    ok = CheckScalar("redMass", redMass, kitty.GetRedMass()) && ok;
    ok = CheckScalar("gridStepX", gridStepX, kitty.GetGridStepX()) && ok;
    ok = CheckScalar("gridStepY", gridStepY, kitty.GetGridStepY()) && ok;
    ok = CheckScalar("gridStepZ", gridStepZ, kitty.GetGridStepZ()) && ok;
    ok = qOutEdges && SameVector(*qOutEdges, MomentumEdgesOut(kitty), "qOutEdges") && ok;
    ok = qSideEdges && SameVector(*qSideEdges, MomentumEdgesSide(kitty), "qSideEdges") && ok;
    ok = qLongEdges && SameVector(*qLongEdges, MomentumEdgesLong(kitty), "qLongEdges") && ok;
    ok = rOutCenters && SameVector(*rOutCenters, GridCentersX(kitty), "rOutCenters") && ok;
    ok = rSideCenters && SameVector(*rSideCenters, GridCentersY(kitty), "rSideCenters") && ok;
    ok = rLongCenters && SameVector(*rLongCenters, GridCentersZ(kitty), "rLongCenters") && ok;
    if(!kernel){
        std::cout << "LoadCATS3DFastKernel: missing kernel branch.\n";
        ok = false;
    }
    else if(kernel->size() != kitty.GetFastKernelStorageSize()){
        std::cout << "CATS3D fast-kernel cache mismatch: kernel size "
                  << kernel->size() << " vs " << kitty.GetFastKernelStorageSize() << "\n";
        ok = false;
    }

    if(!ok){
        std::cout << "LoadCATS3DFastKernel: refused incompatible kernel.\n";
        return false;
    }

    kitty.SetFastWaveFunctionKernel(*kernel, subdivisions);
    std::cout << "Loaded CATS3D fast kernel: " << fileName
              << " entries=" << kernel->size()
              << " qBins=(" << nQOut << "," << nQSide << "," << nQLong << ")"
              << " rBins=(" << nRX << "," << nRY << "," << nRZ << ")"
              << " subdivisions=" << subdivisions << "\n";
    return true;
}

bool SaveCATS3DFastKernelByQ(CATS3D& kitty,
                             const TString& fileName,
                             const TString& treeName,
                             const bool buildIfNeeded)
{
    if(buildIfNeeded && !kitty.IsFastKernelReady()){
        kitty.BuildFastKernelOnly();
    }
    if(!kitty.IsFastKernelReady()){
        std::cout << "SaveCATS3DFastKernelByQ: fast kernel is not ready.\n";
        return false;
    }

    const std::vector<float>& kernel = kitty.GetFastWaveFunctionKernel();
    const size_t expectedSize = kitty.GetFastKernelStorageSize();
    if(kernel.size() != expectedSize){
        std::cout << "SaveCATS3DFastKernelByQ: kernel size mismatch before write.\n";
        return false;
    }

    std::unique_ptr<TFile> fout(TFile::Open(fileName, "RECREATE"));
    if(!fout || fout->IsZombie()){
        std::cout << "SaveCATS3DFastKernelByQ: cannot open " << fileName << "\n";
        return false;
    }

    unsigned version = 1;
    unsigned subdivisions = kitty.GetFastKernelSubdivisions();
    unsigned nQOut = kitty.GetNumMomBinsOut();
    unsigned nQSide = kitty.GetNumMomBinsSide();
    unsigned nQLong = kitty.GetNumMomBinsLong();
    unsigned nRX = kitty.GetNumGridBinsX();
    unsigned nRY = kitty.GetNumGridBinsY();
    unsigned nRZ = kitty.GetNumGridBinsZ();
    unsigned maxPw = kitty.GetMaxPw();
    unsigned wfMode = unsigned(kitty.GetWfMode());
    unsigned identical = kitty.GetIdenticalParticles() ? 1u : 0u;
    int q1q2 = kitty.GetQ1Q2();
    double redMass = kitty.GetRedMass();
    double gridStepX = kitty.GetGridStepX();
    double gridStepY = kitty.GetGridStepY();
    double gridStepZ = kitty.GetGridStepZ();

    std::vector<double> qOutEdges = MomentumEdgesOut(kitty);
    std::vector<double> qSideEdges = MomentumEdgesSide(kitty);
    std::vector<double> qLongEdges = MomentumEdgesLong(kitty);
    std::vector<double> rOutCenters = GridCentersX(kitty);
    std::vector<double> rSideCenters = GridCentersY(kitty);
    std::vector<double> rLongCenters = GridCentersZ(kitty);

    const TString metaName = treeName + "_meta";
    TTree meta(metaName, metaName);
    meta.Branch("version", &version);
    meta.Branch("subdivisions", &subdivisions);
    meta.Branch("nQOut", &nQOut);
    meta.Branch("nQSide", &nQSide);
    meta.Branch("nQLong", &nQLong);
    meta.Branch("nRX", &nRX);
    meta.Branch("nRY", &nRY);
    meta.Branch("nRZ", &nRZ);
    meta.Branch("maxPw", &maxPw);
    meta.Branch("wfMode", &wfMode);
    meta.Branch("identical", &identical);
    meta.Branch("q1q2", &q1q2);
    meta.Branch("redMass", &redMass);
    meta.Branch("gridStepX", &gridStepX);
    meta.Branch("gridStepY", &gridStepY);
    meta.Branch("gridStepZ", &gridStepZ);
    meta.Branch("qOutEdges", &qOutEdges);
    meta.Branch("qSideEdges", &qSideEdges);
    meta.Branch("qLongEdges", &qLongEdges);
    meta.Branch("rOutCenters", &rOutCenters);
    meta.Branch("rSideCenters", &rSideCenters);
    meta.Branch("rLongCenters", &rLongCenters);
    meta.Fill();
    meta.Write("", TObject::kOverwrite);

    const size_t nGrid = kitty.GetFastKernelGridSize();
    unsigned qIndex = 0;
    unsigned iqOut = 0;
    unsigned iqSide = 0;
    unsigned iqLong = 0;
    double qOut = 0.0;
    double qSide = 0.0;
    double qLong = 0.0;
    std::vector<float> kernelR;
    kernelR.reserve(nGrid);

    TTree tree(treeName, treeName);
    tree.Branch("qIndex", &qIndex);
    tree.Branch("iqOut", &iqOut);
    tree.Branch("iqSide", &iqSide);
    tree.Branch("iqLong", &iqLong);
    tree.Branch("qOut", &qOut);
    tree.Branch("qSide", &qSide);
    tree.Branch("qLong", &qLong);
    tree.Branch("kernelR", &kernelR);

    for(iqOut=0; iqOut<nQOut; ++iqOut){
        qOut = kitty.GetMomentumOut(iqOut);
        for(iqSide=0; iqSide<nQSide; ++iqSide){
            qSide = kitty.GetMomentumSide(iqSide);
            for(iqLong=0; iqLong<nQLong; ++iqLong){
                qLong = kitty.GetMomentumLong(iqLong);
                qIndex = (iqOut*nQSide + iqSide)*nQLong + iqLong;
                const size_t offset = size_t(qIndex)*nGrid;
                kernelR.assign(kernel.begin()+offset, kernel.begin()+offset+nGrid);
                tree.Fill();
            }
        }
    }
    tree.Write("", TObject::kOverwrite);
    fout->Close();

    std::cout << "Saved CATS3D fast kernel by q: " << fileName
              << " qEntries=" << size_t(nQOut)*size_t(nQSide)*size_t(nQLong)
              << " rPerEntry=" << nGrid
              << " totalFloats=" << kernel.size()
              << " qBins=(" << nQOut << "," << nQSide << "," << nQLong << ")"
              << " rBins=(" << nRX << "," << nRY << "," << nRZ << ")"
              << " subdivisions=" << subdivisions << "\n";
    return true;
}

bool LoadCATS3DFastKernelByQ(CATS3D& kitty,
                             const TString& fileName,
                             const TString& treeName)
{
    std::unique_ptr<TFile> fin(TFile::Open(fileName, "READ"));
    if(!fin || fin->IsZombie()){
        std::cout << "LoadCATS3DFastKernelByQ: cannot open " << fileName << "\n";
        return false;
    }

    const TString metaName = treeName + "_meta";
    TTree* meta = dynamic_cast<TTree*>(fin->Get(metaName));
    TTree* tree = dynamic_cast<TTree*>(fin->Get(treeName));
    if(!meta || !tree){
        std::cout << "LoadCATS3DFastKernelByQ: missing tree(s) "
                  << metaName << " or " << treeName << " in " << fileName << "\n";
        return false;
    }
    if(meta->GetEntries() <= 0 || tree->GetEntries() <= 0){
        std::cout << "LoadCATS3DFastKernelByQ: empty tree.\n";
        return false;
    }

    unsigned version = 0;
    unsigned subdivisions = 0;
    unsigned nQOut = 0;
    unsigned nQSide = 0;
    unsigned nQLong = 0;
    unsigned nRX = 0;
    unsigned nRY = 0;
    unsigned nRZ = 0;
    unsigned maxPw = 0;
    unsigned wfMode = 0;
    unsigned identical = 0;
    int q1q2 = 0;
    double redMass = 0.0;
    double gridStepX = 0.0;
    double gridStepY = 0.0;
    double gridStepZ = 0.0;
    std::vector<double>* qOutEdges = nullptr;
    std::vector<double>* qSideEdges = nullptr;
    std::vector<double>* qLongEdges = nullptr;
    std::vector<double>* rOutCenters = nullptr;
    std::vector<double>* rSideCenters = nullptr;
    std::vector<double>* rLongCenters = nullptr;

    meta->SetBranchAddress("version", &version);
    meta->SetBranchAddress("subdivisions", &subdivisions);
    meta->SetBranchAddress("nQOut", &nQOut);
    meta->SetBranchAddress("nQSide", &nQSide);
    meta->SetBranchAddress("nQLong", &nQLong);
    meta->SetBranchAddress("nRX", &nRX);
    meta->SetBranchAddress("nRY", &nRY);
    meta->SetBranchAddress("nRZ", &nRZ);
    meta->SetBranchAddress("maxPw", &maxPw);
    meta->SetBranchAddress("wfMode", &wfMode);
    meta->SetBranchAddress("identical", &identical);
    meta->SetBranchAddress("q1q2", &q1q2);
    meta->SetBranchAddress("redMass", &redMass);
    meta->SetBranchAddress("gridStepX", &gridStepX);
    meta->SetBranchAddress("gridStepY", &gridStepY);
    meta->SetBranchAddress("gridStepZ", &gridStepZ);
    meta->SetBranchAddress("qOutEdges", &qOutEdges);
    meta->SetBranchAddress("qSideEdges", &qSideEdges);
    meta->SetBranchAddress("qLongEdges", &qLongEdges);
    meta->SetBranchAddress("rOutCenters", &rOutCenters);
    meta->SetBranchAddress("rSideCenters", &rSideCenters);
    meta->SetBranchAddress("rLongCenters", &rLongCenters);
    meta->GetEntry(0);

    bool ok = true;
    ok = CheckScalar("version", version, 1u) && ok;
    ok = CheckScalar("subdivisions", subdivisions, kitty.GetFastKernelSubdivisions()) && ok;
    ok = CheckScalar("nQOut", nQOut, kitty.GetNumMomBinsOut()) && ok;
    ok = CheckScalar("nQSide", nQSide, kitty.GetNumMomBinsSide()) && ok;
    ok = CheckScalar("nQLong", nQLong, kitty.GetNumMomBinsLong()) && ok;
    ok = CheckScalar("nRX", nRX, kitty.GetNumGridBinsX()) && ok;
    ok = CheckScalar("nRY", nRY, kitty.GetNumGridBinsY()) && ok;
    ok = CheckScalar("nRZ", nRZ, kitty.GetNumGridBinsZ()) && ok;
    ok = CheckScalar("maxPw", maxPw, unsigned(kitty.GetMaxPw())) && ok;
    ok = CheckScalar("wfMode", wfMode, unsigned(kitty.GetWfMode())) && ok;
    ok = CheckScalar("identical", identical, kitty.GetIdenticalParticles() ? 1u : 0u) && ok;
    if(q1q2 != kitty.GetQ1Q2()){
        std::cout << "CATS3D fast-kernel cache mismatch: q1q2 "
                  << q1q2 << " vs " << kitty.GetQ1Q2() << "\n";
        ok = false;
    }
    ok = CheckScalar("redMass", redMass, kitty.GetRedMass()) && ok;
    ok = CheckScalar("gridStepX", gridStepX, kitty.GetGridStepX()) && ok;
    ok = CheckScalar("gridStepY", gridStepY, kitty.GetGridStepY()) && ok;
    ok = CheckScalar("gridStepZ", gridStepZ, kitty.GetGridStepZ()) && ok;
    ok = qOutEdges && SameVector(*qOutEdges, MomentumEdgesOut(kitty), "qOutEdges") && ok;
    ok = qSideEdges && SameVector(*qSideEdges, MomentumEdgesSide(kitty), "qSideEdges") && ok;
    ok = qLongEdges && SameVector(*qLongEdges, MomentumEdgesLong(kitty), "qLongEdges") && ok;
    ok = rOutCenters && SameVector(*rOutCenters, GridCentersX(kitty), "rOutCenters") && ok;
    ok = rSideCenters && SameVector(*rSideCenters, GridCentersY(kitty), "rSideCenters") && ok;
    ok = rLongCenters && SameVector(*rLongCenters, GridCentersZ(kitty), "rLongCenters") && ok;

    const size_t nMom = kitty.GetFastKernelMomSize();
    const size_t nGrid = kitty.GetFastKernelGridSize();
    if(size_t(tree->GetEntries()) != nMom){
        std::cout << "CATS3D fast-kernel cache mismatch: q entries "
                  << tree->GetEntries() << " vs " << nMom << "\n";
        ok = false;
    }
    if(!ok){
        std::cout << "LoadCATS3DFastKernelByQ: refused incompatible kernel.\n";
        return false;
    }

    unsigned qIndex = 0;
    unsigned iqOut = 0;
    unsigned iqSide = 0;
    unsigned iqLong = 0;
    double qOut = 0.0;
    double qSide = 0.0;
    double qLong = 0.0;
    std::vector<float>* kernelR = nullptr;
    tree->SetBranchAddress("qIndex", &qIndex);
    tree->SetBranchAddress("iqOut", &iqOut);
    tree->SetBranchAddress("iqSide", &iqSide);
    tree->SetBranchAddress("iqLong", &iqLong);
    tree->SetBranchAddress("qOut", &qOut);
    tree->SetBranchAddress("qSide", &qSide);
    tree->SetBranchAddress("qLong", &qLong);
    tree->SetBranchAddress("kernelR", &kernelR);

    std::vector<float> kernel(nMom*nGrid, 0.0f);
    for(Long64_t iEntry=0; iEntry<tree->GetEntries(); ++iEntry){
        tree->GetEntry(iEntry);
        const unsigned expectedQIndex = (iqOut*nQSide + iqSide)*nQLong + iqLong;
        if(qIndex != expectedQIndex || qIndex >= nMom){
            std::cout << "LoadCATS3DFastKernelByQ: bad qIndex at entry "
                      << iEntry << "\n";
            return false;
        }
        if(!kernelR || kernelR->size() != nGrid){
            std::cout << "LoadCATS3DFastKernelByQ: bad kernelR size at qIndex "
                      << qIndex << "\n";
            return false;
        }
        std::copy(kernelR->begin(), kernelR->end(),
                  kernel.begin()+size_t(qIndex)*nGrid);
    }

    kitty.SetFastWaveFunctionKernel(kernel, subdivisions);
    std::cout << "Loaded CATS3D fast kernel by q: " << fileName
              << " qEntries=" << tree->GetEntries()
              << " rPerEntry=" << nGrid
              << " totalFloats=" << kernel.size()
              << " qBins=(" << nQOut << "," << nQSide << "," << nQLong << ")"
              << " rBins=(" << nRX << "," << nRY << "," << nRZ << ")"
              << " subdivisions=" << subdivisions << "\n";
    return true;
}
