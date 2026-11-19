#include "TFile.h"
#include "TH3D.h"
#include "TString.h"
#include "TTree.h"

#include <cmath>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <vector>

namespace {

// Same convention as the current CATS3D fast-kernel storage:
// rIndex = (iRout*nRside + iRside)*nRlong + iRlong.
size_t RIndex(const unsigned iRout,
              const unsigned iRside,
              const unsigned iRlong,
              const unsigned nRside,
              const unsigned nRlong)
{
    return (size_t(iRout)*nRside + iRside)*nRlong + iRlong;
}

// Same convention as the current CATS3D fast-kernel storage:
// qIndex = (iQout*nQside + iQside)*nQlong + iQlong.
unsigned QIndex(const unsigned iQout,
                const unsigned iQside,
                const unsigned iQlong,
                const unsigned nQside,
                const unsigned nQlong)
{
    return (iQout*nQside + iQside)*nQlong + iQlong;
}

// 3D Gaussian source in the same "pair relative source" convention used for
// the CATS3D convolution. The absolute normalization is not important if we
// divide by sourceInt after the integral.
double GaussianSource3D(const double rOut,
                        const double rSide,
                        const double rLong,
                        const double Rout,
                        const double Rside,
                        const double Rlong)
{
    if(Rout <= 0.0 || Rside <= 0.0 || Rlong <= 0.0) return 0.0;

    const double x = rOut/Rout;
    const double y = rSide/Rside;
    const double z = rLong/Rlong;
    return std::exp(-0.25*(x*x + y*y + z*z));
}

} // namespace

void read_wavefunction_kernel_example(
    const char* kernelFileName = "./OutputFiles/kernel_cache_coulomb.root",
    const char* treeNameInput = "CATS3DFastKernel6DByQ")
{
    const TString treeName(treeNameInput);
    const TString metaName = treeName + "_meta";

    std::unique_ptr<TFile> fin(TFile::Open(kernelFileName, "READ"));
    if(!fin || fin->IsZombie()){
        throw std::runtime_error(Form("cannot open kernel file: %s", kernelFileName));
    }

    // The by-q file has two TTrees:
    //   1. CATS3DFastKernel6DByQ_meta:
    //      one entry containing grid sizes, q-bin edges, r-bin centers,
    //      CATS settings, and r-grid step sizes.
    //   2. CATS3DFastKernel6DByQ:
    //      one entry per q-bin. Each entry contains one vector<float> kernelR,
    //      i.e. |Psi(q,r)|^2 evaluated on all 3D r-grid points.
    TTree* meta = dynamic_cast<TTree*>(fin->Get(metaName));
    TTree* tree = dynamic_cast<TTree*>(fin->Get(treeName));
    if(!meta || !tree){
        throw std::runtime_error(Form("missing %s or %s in %s",
                                      metaName.Data(), treeName.Data(), kernelFileName));
    }

    // Metadata branches. q edges are in MeV/c, r centers are in fm.
    unsigned nQOut = 0;
    unsigned nQSide = 0;
    unsigned nQLong = 0;
    unsigned nRout = 0;
    unsigned nRside = 0;
    unsigned nRlong = 0;
    double gridStepOut = 0.0;
    double gridStepSide = 0.0;
    double gridStepLong = 0.0;
    std::vector<double>* qOutEdges = nullptr;
    std::vector<double>* qSideEdges = nullptr;
    std::vector<double>* qLongEdges = nullptr;
    std::vector<double>* rOutCenters = nullptr;
    std::vector<double>* rSideCenters = nullptr;
    std::vector<double>* rLongCenters = nullptr;

    meta->SetBranchAddress("nQOut", &nQOut);
    meta->SetBranchAddress("nQSide", &nQSide);
    meta->SetBranchAddress("nQLong", &nQLong);
    meta->SetBranchAddress("nRX", &nRout);
    meta->SetBranchAddress("nRY", &nRside);
    meta->SetBranchAddress("nRZ", &nRlong);
    meta->SetBranchAddress("gridStepX", &gridStepOut);
    meta->SetBranchAddress("gridStepY", &gridStepSide);
    meta->SetBranchAddress("gridStepZ", &gridStepLong);
    meta->SetBranchAddress("qOutEdges", &qOutEdges);
    meta->SetBranchAddress("qSideEdges", &qSideEdges);
    meta->SetBranchAddress("qLongEdges", &qLongEdges);
    meta->SetBranchAddress("rOutCenters", &rOutCenters);
    meta->SetBranchAddress("rSideCenters", &rSideCenters);
    meta->SetBranchAddress("rLongCenters", &rLongCenters);
    meta->GetEntry(0);

    if(!qOutEdges || !qSideEdges || !qLongEdges ||
       !rOutCenters || !rSideCenters || !rLongCenters){
        throw std::runtime_error("metadata vectors are missing");
    }

    const size_t nRGrid = size_t(nRout)*size_t(nRside)*size_t(nRlong);
    const size_t nQGrid = size_t(nQOut)*size_t(nQSide)*size_t(nQLong);

    std::cout << "Loaded kernel metadata\n"
              << "  q bins = (" << nQOut << ", " << nQSide << ", " << nQLong
              << "), total q entries = " << nQGrid << "\n"
              << "  r bins = (" << nRout << ", " << nRside << ", " << nRlong
              << "), r points per q = " << nRGrid << "\n"
              << "  r steps = (" << gridStepOut << ", "
              << gridStepSide << ", " << gridStepLong << ") fm\n";

    // Kernel tree branches. qIndex is the flattened q-bin index. kernelR stores
    // all |Psi(q,r)|^2 values for this q-bin.
    unsigned qIndex = 0;
    unsigned iQout = 0;
    unsigned iQside = 0;
    unsigned iQlong = 0;
    double qOut = 0.0;
    double qSide = 0.0;
    double qLong = 0.0;
    std::vector<float>* kernelR = nullptr;

    tree->SetBranchAddress("qIndex", &qIndex);
    tree->SetBranchAddress("iqOut", &iQout);
    tree->SetBranchAddress("iqSide", &iQside);
    tree->SetBranchAddress("iqLong", &iQlong);
    tree->SetBranchAddress("qOut", &qOut);
    tree->SetBranchAddress("qSide", &qSide);
    tree->SetBranchAddress("qLong", &qLong);
    tree->SetBranchAddress("kernelR", &kernelR);

    // Example 1: read one q-bin and one r-bin.
    // Here we choose the central q-bin and central r-bin.
    const unsigned targetIQout = nQOut/2;
    const unsigned targetIQside = nQSide/2;
    const unsigned targetIQlong = nQLong/2;
    const unsigned targetQIndex =
        QIndex(targetIQout, targetIQside, targetIQlong, nQSide, nQLong);

    tree->GetEntry(targetQIndex);
    if(!kernelR || kernelR->size() != nRGrid){
        throw std::runtime_error("kernelR has wrong size at target q-bin");
    }

    const unsigned targetIRout = nRout/2;
    const unsigned targetIRside = nRside/2;
    const unsigned targetIRlong = nRlong/2;
    const size_t targetRIndex =
        RIndex(targetIRout, targetIRside, targetIRlong, nRside, nRlong);

    std::cout << "\nExample 1: one |Psi|^2 value\n"
              << "  q bin index = " << qIndex
              << ", q = (" << qOut << ", " << qSide << ", " << qLong
              << ") MeV/c\n"
              << "  r bin index = " << targetRIndex
              << ", r = (" << (*rOutCenters)[targetIRout] << ", "
              << (*rSideCenters)[targetIRside] << ", "
              << (*rLongCenters)[targetIRlong] << ") fm\n"
              << "  |Psi(q,r)|^2 = " << (*kernelR)[targetRIndex] << "\n";

    // Example 2: convolve this stored |Psi|^2 kernel with a Gaussian source.
    // This is the essential operation used during the fast CATS3D fit:
    //
    //   C(q) = integral d^3r S(r; Rout,Rside,Rlong) |Psi(q,r)|^2
    //          -----------------------------------------------------
    //          integral d^3r S(r; Rout,Rside,Rlong)
    //
    // Since the r-grid is discrete, d^3r = gridStepOut*gridStepSide*gridStepLong.
    const double Rout = 5.0;
    const double Rside = 5.0;
    const double Rlong = 5.0;
    const double dV = gridStepOut*gridStepSide*gridStepLong;

    std::vector<double> sourceWeight(nRGrid, 0.0);
    double sourceInt = 0.0;
    for(unsigned iRout=0; iRout<nRout; ++iRout){
        for(unsigned iRside=0; iRside<nRside; ++iRside){
            for(unsigned iRlong=0; iRlong<nRlong; ++iRlong){
                const size_t rIndex = RIndex(iRout, iRside, iRlong,
                                             nRside, nRlong);
                const double src = GaussianSource3D((*rOutCenters)[iRout],
                                                    (*rSideCenters)[iRside],
                                                    (*rLongCenters)[iRlong],
                                                    Rout, Rside, Rlong);
                sourceWeight[rIndex] = src*dV;
                sourceInt += sourceWeight[rIndex];
            }
        }
    }

    // Build a TH3D C(qout,qside,qlong) using the same q-bin edges as the kernel.
    auto hCk = std::make_unique<TH3D>(
        "hCk_fromKernel_example",
        "C(q) from cached |#Psi|^{2};q_{out} (MeV/c);q_{side} (MeV/c);q_{long} (MeV/c)",
        nQOut, qOutEdges->data(),
        nQSide, qSideEdges->data(),
        nQLong, qLongEdges->data());
    hCk->SetDirectory(nullptr);

    for(Long64_t iEntry=0; iEntry<tree->GetEntries(); ++iEntry){
        tree->GetEntry(iEntry);
        if(!kernelR || kernelR->size() != nRGrid){
            throw std::runtime_error(Form("bad kernelR at entry %lld", iEntry));
        }

        double numerator = 0.0;
        for(size_t iR=0; iR<nRGrid; ++iR){
            numerator += sourceWeight[iR]*double((*kernelR)[iR]);
        }
        const double ck = numerator/sourceInt;

        // ROOT histogram bins are 1-based, while C++ indices here are 0-based.
        hCk->SetBinContent(iQout+1, iQside+1, iQlong+1, ck);
        hCk->SetBinError(iQout+1, iQside+1, iQlong+1, 0.0);
    }

    TFile fout("read_wavefunction_kernel_example_output.root", "RECREATE");
    hCk->Write();
    fout.Close();

    std::cout << "\nExample 2: Gaussian convolution finished\n"
              << "  R = (" << Rout << ", " << Rside << ", " << Rlong << ") fm\n"
              << "  source integral on finite r-grid = " << sourceInt << "\n"
              << "  wrote read_wavefunction_kernel_example_output.root\n";
}
