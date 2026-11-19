#include "Basics3D.h"
#include "/Users/wuwenya/Documents/code/femto_flow/CATS/CATS_install/install/include/CATSconstants.h"

#include "TCanvas.h"
#include "TFile.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TStyle.h"

#include <array>
#include <iomanip>
#include <iostream>
#include <memory>

namespace {
constexpr double lambdaGen = 0.83;
constexpr double weakFraction = 0.14;
constexpr double lambdaEffective = lambdaGen*(1.0-weakFraction)*(1.0-weakFraction);
constexpr std::array<double,3> radii = {8.64286, 7.52793, 9.02600};

void Configure(CATS3D& kitty, unsigned nOut, double outMin, double outMax,
               unsigned nSide, double sideMin, double sideMax,
               unsigned nLong, double longMin, double longMax) {
    kitty.SetMomBinsOut(nOut, outMin, outMax);
    kitty.SetMomBinsSide(nSide, sideMin, sideMax);
    kitty.SetMomBinsLong(nLong, longMin, longMax);
    kitty.SetGridX(31, -50.0, 50.0);
    kitty.SetGridY(31, -50.0, 50.0);
    kitty.SetGridZ(31, -50.0, 50.0);
    kitty.SetFastKernelMode(true);
    kitty.SetFastKernelSubdivisions(1);
    CATSparameters pars(CATSparameters::tSource3D, 3, true);
    for(unsigned i=0;i<3;++i) pars.SetParameter(i,radii[i]);
    kitty.SetAnaSource(Basics_Source_Gauss3D,pars);
    kitty.SetPdgId(211,211);
    kitty.SetQuantumStatistics(1);
    kitty.SetQ1Q2(1);
    kitty.SetRedMass(0.5*Mass_Pic);
    kitty.SetWfMode(CATS3D::wfCoulomb);
}

double Total(const CATS3D& kitty, unsigned out, unsigned side, unsigned lon) {
    return 1.0+lambdaEffective*(kitty.GetCorrFun(out,side,lon)-1.0);
}
}

int main() {
    gStyle->SetOptStat(0);
    std::cout << std::setprecision(8)
              << "lambdaGen=" << lambdaGen << " weak=" << weakFraction
              << " lambdaEffective=" << lambdaEffective << " R="
              << radii[0] << ',' << radii[1] << ',' << radii[2] << " fm\n";

    CATS3D coarseAxis;
    Configure(coarseAxis,20,-200,200,1,-10,10,1,-10,10);
    coarseAxis.KillTheCat();
    CATS3D fineAxis;
    Configure(fineAxis,40,-200,200,1,-10,10,1,-10,10);
    fineAxis.KillTheCat();

    TH1D h20("hCATS_q20_center",";q_{out} (GeV/c);C_{3D}(q_{out},0,0)",20,-.2,.2);
    TH1D h10("hCATS_q10_center",";q_{out} (GeV/c);C_{3D}(q_{out},0,0)",40,-.2,.2);
    TH1D h10avg("hCATS_q10_avg_to_q20",";q_{out} (GeV/c);C_{3D}(q_{out},0,0)",20,-.2,.2);
    for(unsigned i=0;i<20;++i) h20.SetBinContent(i+1,Total(coarseAxis,i,0,0));
    for(unsigned i=0;i<40;++i) h10.SetBinContent(i+1,Total(fineAxis,i,0,0));
    for(unsigned i=0;i<20;++i) h10avg.SetBinContent(i+1,
        0.5*(h10.GetBinContent(2*i+1)+h10.GetBinContent(2*i+2)));

    // Match one positive central 20-MeV 3D voxel [0,20]^3 with the
    // eight corresponding 10-MeV voxel centers.
    CATS3D coarseCube;
    Configure(coarseCube,2,0,40,2,0,40,2,0,40);
    coarseCube.KillTheCat();
    CATS3D fineCube;
    Configure(fineCube,4,0,40,4,0,40,4,0,40);
    fineCube.KillTheCat();
    const double coarseCenter=Total(coarseCube,0,0,0);
    double fineAverage=0;
    for(unsigned ix=0;ix<2;++ix)
      for(unsigned iy=0;iy<2;++iy)
        for(unsigned iz=0;iz<2;++iz)
          fineAverage+=Total(fineCube,ix,iy,iz)/8.0;
    CATS3D finerCube;
    Configure(finerCube,8,0,40,8,0,40,8,0,40);
    finerCube.KillTheCat();
    double finerAverage=0;
    for(unsigned ix=0;ix<4;++ix)
      for(unsigned iy=0;iy<4;++iy)
        for(unsigned iz=0;iz<4;++iz)
          finerAverage+=Total(finerCube,ix,iy,iz)/64.0;
    std::cout << "central 3D voxel [0,20]^3 MeV/c: 20-MeV center="
              << coarseCenter << " average of eight 10-MeV centers="
              << fineAverage << " average of sixty-four 5-MeV centers="
              << finerAverage << " delta(5-10)=" << finerAverage-fineAverage << '\n';
    std::cout << "axis [0,20] MeV/c at qside=qlong=0: 20-MeV center="
              << h20.GetBinContent(11) << " mean of 10-MeV centers="
              << h10avg.GetBinContent(11) << '\n';

    TFile file("derivedFull_Bminus_cent0to10_anyCPR_lambda083_test/cats3d_qbin_comparison.root","RECREATE");
    h20.Write(); h10.Write(); h10avg.Write();
    TH1D summary("hCentralVoxelComparison","[0,20]^{3} MeV/c central voxel;case;C",3,0,3);
    summary.GetXaxis()->SetBinLabel(1,"20-MeV center");
    summary.GetXaxis()->SetBinLabel(2,"10-MeV average");
    summary.GetXaxis()->SetBinLabel(3,"5-MeV average");
    summary.SetBinContent(1,coarseCenter); summary.SetBinContent(2,fineAverage);
    summary.SetBinContent(3,finerAverage);
    summary.Write();
    TCanvas canvas("cMomentumBinComparison","",1000,650);
    canvas.SetLeftMargin(.12); canvas.SetRightMargin(.04); canvas.SetBottomMargin(.13);
    h20.SetLineColor(kBlue+1); h20.SetLineWidth(2);
    h10.SetLineColor(kGray+2); h10.SetLineWidth(2);
    h10avg.SetLineColor(kRed+1); h10avg.SetLineWidth(3);
    h20.GetXaxis()->SetRangeUser(-.12,.12);
    h20.SetMinimum(.96); h20.SetMaximum(1.75);
    h20.Draw("HIST"); h10.Draw("HIST SAME"); h10avg.Draw("HIST SAME");
    TLegend legend(.53,.69,.94,.88);
    legend.SetBorderSize(0); legend.SetFillStyle(0);
    legend.AddEntry(&h20,"20 MeV/c: bin center","l");
    legend.AddEntry(&h10,"10 MeV/c: bin center","l");
    legend.AddEntry(&h10avg,"10 MeV/c averaged to 20 MeV/c","l");
    legend.Draw();
    TLatex label; label.SetNDC(); label.SetTextSize(.035);
    label.DrawLatex(.14,.91,"CATS3D, Coulomb + Gaussian source; fixed R and #lambda");
    canvas.Write();
    canvas.SaveAs("derivedFull_Bminus_cent0to10_anyCPR_lambda083_test/cats3d_qbin_comparison.png");
}
