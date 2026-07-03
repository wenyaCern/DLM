#include "TFile.h"
#include "TH3D.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TString.h"
#include "TStyle.h"
#include "TAxis.h"
#include "TLine.h"
#include <iostream>

void project_CK3D(const char* infile="test_CATS3D_wtCoul_highl.root",
                  const char* hname="hCk3D_plane",
                  double qCut=5.0) {

    gStyle->SetOptStat(0);

    TFile* f = TFile::Open(infile, "READ");
    if(!f || f->IsZombie()){
        std::cout << "Cannot open file: " << infile << std::endl;
        return;
    }

    TH3D* h3 = (TH3D*)f->Get(hname);
    if(!h3){
        std::cout << "Cannot find histogram: " << hname << std::endl;
        f->Close();
        return;
    }

    // -------------------------
    // 1) full projections
    // -------------------------
    TH1D* hQout_full  = (TH1D*)h3->ProjectionX("hQout_full");
    TH1D* hQside_full = (TH1D*)h3->ProjectionY("hQside_full");
    TH1D* hQlong_full = (TH1D*)h3->ProjectionZ("hQlong_full");

    hQout_full ->SetTitle("Full projection on q_{out};q_{out} (MeV/c);Counts");
    hQside_full->SetTitle("Full projection on q_{side};q_{side} (MeV/c);Counts");
    hQlong_full->SetTitle("Full projection on q_{long};q_{long} (MeV/c);Counts");

    // -------------------------
    // 2) central-window projections
    //    |other two coordinates| < qCut
    // -------------------------
    int yBinMin = h3->GetYaxis()->FindBin(-qCut + 1e-6);
    int yBinMax = h3->GetYaxis()->FindBin( qCut - 1e-6);
    int zBinMin = h3->GetZaxis()->FindBin(-qCut + 1e-6);
    int zBinMax = h3->GetZaxis()->FindBin( qCut - 1e-6);
    int xBinMin = h3->GetXaxis()->FindBin(-qCut + 1e-6);
    int xBinMax = h3->GetXaxis()->FindBin( qCut - 1e-6);

    int outNumBin =  (yBinMax - yBinMin + 1)*(zBinMax - zBinMin + 1);
    int sideNumBin = (xBinMax - xBinMin + 1)*(zBinMax - zBinMin + 1);
    int longNumBin = (xBinMax - xBinMin + 1)*(yBinMax - yBinMin + 1);

    TH1D* hQout_cut  = (TH1D*)h3->ProjectionX("hQout_cut",  yBinMin, yBinMax, zBinMin, zBinMax);
    TH1D* hQside_cut = (TH1D*)h3->ProjectionY("hQside_cut", xBinMin, xBinMax, zBinMin, zBinMax);
    TH1D* hQlong_cut = (TH1D*)h3->ProjectionZ("hQlong_cut", xBinMin, xBinMax, yBinMin, yBinMax);

    hQout_cut ->SetTitle(Form("q_{out} projection, |q_{side}|,|q_{long}|<%.1f MeV/c;q_{out} (MeV/c);Counts", qCut));
    hQside_cut->SetTitle(Form("q_{side} projection, |q_{out}|,|q_{long}|<%.1f MeV/c;q_{side} (MeV/c);Counts", qCut));
    hQlong_cut->SetTitle(Form("q_{long} projection, |q_{out}|,|q_{side}|<%.1f MeV/c;q_{long} (MeV/c);Counts", qCut));

    hQout_cut->Scale(1./outNumBin);
    hQside_cut->Scale(1./sideNumBin);
    hQlong_cut->Scale(1./longNumBin);

    // -------------------------
    // draw
    // -------------------------
    TCanvas* c1 = new TCanvas("cProjFull","Full projections",1200,400);
    c1->Divide(3,1);
    c1->cd(1); hQout_full->Draw("hist");
    c1->cd(2); hQside_full->Draw("hist");
    c1->cd(3); hQlong_full->Draw("hist");

    TCanvas* c2 = new TCanvas("cProjCut","Central-window projections",1200,400);
    c2->Divide(3,1);
    c2->cd(1); hQout_cut->Draw("hist");
    c2->cd(2); hQside_cut->Draw("hist");
    c2->cd(3); hQlong_cut->Draw("hist");

    // optional: save
    TFile* fout = TFile::Open("proj_CK3D_wtCoul_highl.root","RECREATE");
    hQout_full->Write();
    hQside_full->Write();
    hQlong_full->Write();
    hQout_cut->Write();
    hQside_cut->Write();
    hQlong_cut->Write();
    fout->Close();

    f->Close();
}