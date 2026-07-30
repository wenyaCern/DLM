#include "Basics3D.h"

#include "TFile.h"
#include "TH1D.h"
#include "TH3D.h"
#include "TString.h"
#include "TSystem.h"
#include "TTree.h"

#include <cstdlib>
#include <cmath>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <vector>

namespace {

double GaussSource3D(const double rOut,
                     const double rSide,
                     const double rLong,
                     const double Rout,
                     const double Rside,
                     const double Rlong)
{
    return 1.0/(std::pow(4.0*Pi, 1.5)*Rout*Rside*Rlong) *
           std::exp(-rOut*rOut/(4.0*Rout*Rout)
                    -rSide*rSide/(4.0*Rside*Rside)
                    -rLong*rLong/(4.0*Rlong*Rlong));
}

TH1D* ProjectCkFemtoSlice1D(TH3D* h3,
                            const TString& axis,
                            const TString& name,
                            const double cutAbsMeV)
{
    if(!h3) return nullptr;
    auto centeredBinRange = [](TAxis* axisIn, const double cutAbs,
                               int& firstBin, int& lastBin) {
        firstBin = 1;
        lastBin = 0;
        for(int iBin=1; iBin<=axisIn->GetNbins(); ++iBin){
            if(std::fabs(axisIn->GetBinCenter(iBin)) <= cutAbs + 1e-12){
                if(lastBin == 0) firstBin = iBin;
                lastBin = iBin;
            }
        }
        if(lastBin == 0){
            const int centerBin = axisIn->FindBin(0.0);
            firstBin = centerBin;
            lastBin = centerBin;
        }
    };

    TH1D* h1 = nullptr;
    if(axis=="out"){
        int y1, y2, z1, z2;
        centeredBinRange(h3->GetYaxis(), cutAbsMeV, y1, y2);
        centeredBinRange(h3->GetZaxis(), cutAbsMeV, z1, z2);
        h1 = dynamic_cast<TH1D*>(h3->ProjectionX(name, y1, y2, z1, z2, "e"));
        h1->Scale(1./double((y2-y1+1)*(z2-z1+1)));
    }
    else if(axis=="side"){
        int x1, x2, z1, z2;
        centeredBinRange(h3->GetXaxis(), cutAbsMeV, x1, x2);
        centeredBinRange(h3->GetZaxis(), cutAbsMeV, z1, z2);
        h1 = dynamic_cast<TH1D*>(h3->ProjectionY(name, x1, x2, z1, z2, "e"));
        h1->Scale(1./double((x2-x1+1)*(z2-z1+1)));
    }
    else if(axis=="long"){
        int x1, x2, y1, y2;
        centeredBinRange(h3->GetXaxis(), cutAbsMeV, x1, x2);
        centeredBinRange(h3->GetYaxis(), cutAbsMeV, y1, y2);
        h1 = dynamic_cast<TH1D*>(h3->ProjectionZ(name, x1, x2, y1, y2, "e"));
        h1->Scale(1./double((x2-x1+1)*(y2-y1+1)));
    }
    if(h1) h1->SetDirectory(nullptr);
    return h1;
}

} // namespace

int main(int argc, char** argv)
{
    TString kernelFile = "../OutputFiles/dataKernel_qData_r31_coulomb_byQ.root";
    if(argc > 1) kernelFile = argv[1];

    double Rout = 5.0;
    double Rside = 5.0;
    double Rlong = 5.0;
    if(argc > 2) Rout = std::atof(argv[2]);
    if(argc > 3) Rside = std::atof(argv[3]);
    if(argc > 4) Rlong = std::atof(argv[4]);

    const TString treeName = "CATS3DFastKernel6DByQ";
    const TString metaName = treeName + "_meta";

    std::unique_ptr<TFile> fin(TFile::Open(kernelFile, "READ"));
    if(!fin || fin->IsZombie()){
        throw std::runtime_error(Form("cannot open kernel file '%s'", kernelFile.Data()));
    }

    TTree* meta = dynamic_cast<TTree*>(fin->Get(metaName));
    TTree* tree = dynamic_cast<TTree*>(fin->Get(treeName));
    if(!meta || !tree){
        throw std::runtime_error("missing by-q kernel tree or metadata tree");
    }

    unsigned nQOut = 0;
    unsigned nQSide = 0;
    unsigned nQLong = 0;
    unsigned nRX = 0;
    unsigned nRY = 0;
    unsigned nRZ = 0;
    double gridStepX = 0.0;
    double gridStepY = 0.0;
    double gridStepZ = 0.0;
    std::vector<double>* qOutEdges = nullptr;
    std::vector<double>* qSideEdges = nullptr;
    std::vector<double>* qLongEdges = nullptr;
    std::vector<double>* rOutCenters = nullptr;
    std::vector<double>* rSideCenters = nullptr;
    std::vector<double>* rLongCenters = nullptr;

    meta->SetBranchAddress("nQOut", &nQOut);
    meta->SetBranchAddress("nQSide", &nQSide);
    meta->SetBranchAddress("nQLong", &nQLong);
    meta->SetBranchAddress("nRX", &nRX);
    meta->SetBranchAddress("nRY", &nRY);
    meta->SetBranchAddress("nRZ", &nRZ);
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

    if(!qOutEdges || !qSideEdges || !qLongEdges ||
       !rOutCenters || !rSideCenters || !rLongCenters){
        throw std::runtime_error("missing metadata vectors");
    }

    const size_t nGrid = size_t(nRX)*size_t(nRY)*size_t(nRZ);
    std::vector<double> sourceWeight(nGrid, 0.0);
    const double cellVolume = gridStepX*gridStepY*gridStepZ;
    double sourceInt = 0.0;
    for(unsigned ix=0; ix<nRX; ++ix){
        for(unsigned iy=0; iy<nRY; ++iy){
            for(unsigned iz=0; iz<nRZ; ++iz){
                const size_t rIndex = (size_t(ix)*nRY + iy)*nRZ + iz;
                const double src = GaussSource3D((*rOutCenters)[ix],
                                                 (*rSideCenters)[iy],
                                                 (*rLongCenters)[iz],
                                                 Rout, Rside, Rlong);
                const double weight = src*cellVolume;
                sourceWeight[rIndex] = weight;
                sourceInt += weight;
            }
        }
    }

    std::unique_ptr<TH3D> hCk(new TH3D("hCk_fromKernel_gauss_R5",
                                      "C(q) from cached |#Psi|^{2} #otimes Gaussian R=5;q_{out} (MeV/c);q_{side} (MeV/c);q_{long} (MeV/c)",
                                      nQOut, qOutEdges->data(),
                                      nQSide, qSideEdges->data(),
                                      nQLong, qLongEdges->data()));
    hCk->SetDirectory(nullptr);

    unsigned qIndex = 0;
    unsigned iqOut = 0;
    unsigned iqSide = 0;
    unsigned iqLong = 0;
    std::vector<float>* kernelR = nullptr;
    tree->SetBranchAddress("qIndex", &qIndex);
    tree->SetBranchAddress("iqOut", &iqOut);
    tree->SetBranchAddress("iqSide", &iqSide);
    tree->SetBranchAddress("iqLong", &iqLong);
    tree->SetBranchAddress("kernelR", &kernelR);

    for(Long64_t iEntry=0; iEntry<tree->GetEntries(); ++iEntry){
        tree->GetEntry(iEntry);
        if(!kernelR || kernelR->size() != nGrid){
            throw std::runtime_error(Form("bad kernelR at entry %lld", iEntry));
        }
        double corr = 0.0;
        for(size_t ir=0; ir<nGrid; ++ir){
            corr += sourceWeight[ir]*double((*kernelR)[ir]);
        }
        hCk->SetBinContent(iqOut+1, iqSide+1, iqLong+1, corr/sourceInt);
        hCk->SetBinError(iqOut+1, iqSide+1, iqLong+1, 0.0);
    }

    const TString outFile = Form("../OutputFiles/convolved_gauss_R%.1f_%.1f_%.1f_from_%s",
                                 Rout, Rside, Rlong,
                                 gSystem->BaseName(kernelFile));
    std::unique_ptr<TH1D> hOut(ProjectCkFemtoSlice1D(
        hCk.get(), "out", "hCk_fromKernel_gauss_R5_proj_out", 40.0));
    std::unique_ptr<TH1D> hSide(ProjectCkFemtoSlice1D(
        hCk.get(), "side", "hCk_fromKernel_gauss_R5_proj_side", 40.0));
    std::unique_ptr<TH1D> hLong(ProjectCkFemtoSlice1D(
        hCk.get(), "long", "hCk_fromKernel_gauss_R5_proj_long", 40.0));

    TFile fout(outFile, "RECREATE");
    hCk->Write();
    if(hOut) hOut->Write();
    if(hSide) hSide->Write();
    if(hLong) hLong->Write();
    fout.Close();

    std::cout << "Convolved cached |Psi|^2 with Gaussian source:"
              << " R=(" << Rout << "," << Rside << "," << Rlong << ") fm"
              << " sourceInt=" << sourceInt
              << " qEntries=" << tree->GetEntries()
              << " rPerEntry=" << nGrid
              << "\nWrote " << outFile << "\n";
    return 0;
}
