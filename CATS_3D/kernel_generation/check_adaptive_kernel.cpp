#include "AdaptivePipiKernel.h"
#include "TFile.h"
#include "TTree.h"

#include <cmath>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
void CheckMetadata(const char* path)
{
    std::unique_ptr<TFile> file(TFile::Open(path, "READ"));
    if (!file || file->IsZombie()) throw std::runtime_error("Cannot open kernel");
    auto* meta = dynamic_cast<TTree*>(file->Get("KernelMetadata"));
    auto* tree = dynamic_cast<TTree*>(file->Get("KernelByQ"));
    if (!meta || !tree || meta->GetEntries() != 1)
        throw std::runtime_error("Missing kernel metadata or data tree");

    std::vector<double>* qEdges = nullptr;
    std::vector<double>* rEdges = nullptr;
    meta->SetBranchAddress("qOutEdgesMeV", &qEdges);
    meta->SetBranchAddress("rOutEdgesFm", &rEdges);
    meta->GetEntry(0);
    if (!qEdges || !rEdges || qEdges->size() < 2 || rEdges->size() < 2)
        throw std::runtime_error("Missing kernel axes");
    const auto nQ = static_cast<unsigned>(qEdges->size() - 1);
    const auto nR = static_cast<unsigned>(rEdges->size() - 1);
    const auto expected = static_cast<Long64_t>(nQ) * nQ * nQ;
    if (tree->GetEntries() != expected)
        throw std::runtime_error("Incomplete kernel: q-cell count differs from metadata");

    unsigned ix = 0, iy = 0, iz = 0;
    tree->SetBranchStatus("*", 0);
    for (const char* name : {"iqOut", "iqSide", "iqLong"}) tree->SetBranchStatus(name, 1);
    tree->SetBranchAddress("iqOut", &ix);
    tree->SetBranchAddress("iqSide", &iy);
    tree->SetBranchAddress("iqLong", &iz);
    for (Long64_t i = 0; i < expected; ++i) {
        tree->GetEntry(i);
        if (ix != static_cast<unsigned>(i / (nQ*nQ)) ||
            iy != static_cast<unsigned>((i / nQ) % nQ) ||
            iz != static_cast<unsigned>(i % nQ))
            throw std::runtime_error("Missing or out-of-order q cell");
    }

    std::vector<float>* values = nullptr;
    tree->SetBranchStatus("wavefunction2", 1);
    tree->SetBranchAddress("wavefunction2", &values);
    for (const Long64_t i : {Long64_t(0), expected / 2, expected - 1}) {
        tree->GetEntry(i);
        if (!values || values->size() != static_cast<size_t>(nR)*nR*nR)
            throw std::runtime_error("Sampled q cell has an invalid spatial payload");
        for (size_t j : {size_t(0), values->size()/2, values->size()-1})
            if (!std::isfinite((*values)[j]) || (*values)[j] < 0.0f)
                throw std::runtime_error("Sampled wavefunction value is invalid");
    }
    std::cout << "Kernel metadata and q-cell order valid: " << nQ << "^3 x "
              << nR << "^3, " << expected << " q cells\n";
}
}

int main(int argc, char** argv)
{
    if(argc!=2 && argc!=3) return 2;
    try{
        if (argc == 3) {
            if (std::string(argv[2]) != "--metadata") return 2;
            CheckMetadata(argv[1]);
            return 0;
        }
        AdaptivePipiKernel kernel(argv[1]);
        const auto c=kernel.Fold(8.0,8.0,8.0);
        if(c.size()!=1 || !std::isfinite(c[0]) || c[0]<0.9 || c[0]>1.1)
            throw std::runtime_error("Unexpected folded correlation");
        const auto both=kernel.AxisOverlaps(0,-5,5);
        const auto positive=kernel.AxisOverlaps(0,0,5);
        if(both.size()!=1 || positive.size()!=1 ||
           std::abs(both[0].second-10)>1e-9 || std::abs(positive[0].second-5)>1e-9)
            throw std::runtime_error("Positive/negative q bin overlap failed");
        std::cout << "Folded C(q=2.5,2.5,2.5; R=8 fm)=" << c[0]
                  << "; symmetric q-bin overlap passed\n";
        return 0;
    }catch(const std::exception& e){
        std::cerr << e.what() << '\n';
        return 1;
    }
}
