#include "/Users/wuwenya/Documents/code/femto_flow/CATS/CATS_install/install/include/DLM_Source.h"
#include "/Users/wuwenya/Documents/code/femto_flow/CATS/CATS_install/install/include/DLM_Random.h"

#include "TCanvas.h"
#include "TFile.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TMath.h"
#include "TNtuple.h"
#include "TString.h"
#include "TStyle.h"

#include <cmath>
#include <iostream>

namespace {

const TString kBaseDir =
    "/Users/wuwenya/Documents/code/femto_flow/CATS/CATS_install/CATS_Pion_PbPb";
const TString kPRFile = kBaseDir + "/Files/EPOSPbPb_pi_piReso_noOmega_PRF.root";
const TString kRRFile = kBaseDir + "/Files/EPOSPbPb_piReso_piReso_noOmega_PRF.root";
const TString kOutDir = kBaseDir + "/CMake/out/Aug";

double EvalPureGaussSource(const double r, const double radius)
{
  double pars[4] = {0.0, r, 0.0, radius};
  return GaussSource(pars);
}

DLM_CleverMcLevyResoTM* BuildResonanceSource(const double radius,
                                             const double alpha,
                                             const double fracReso,
                                             const double tauReso,
                                             const double massReso,
                                             const double cutoffMeV)
{
  (void)massReso;
  DLM_CleverMcLevyResoTM* source = new DLM_CleverMcLevyResoTM();
  source->InitStability(1, alpha-1e-6, alpha+1e-6);
  source->InitScale(160, 4.0, 30.0);
  source->InitRad(700, 1e-4, 160.0);
  source->InitType(2);
  source->SetUpReso(0, fracReso);
  source->SetUpReso(1, fracReso);

  Float_t k_D = 0;
  Float_t fP1 = 0;
  Float_t fP2 = 0;
  Float_t fM1 = 0;
  Float_t fM2 = 0;
  Float_t Tau1 = 0;
  Float_t Tau2 = 0;
  Float_t AngleRcP1 = 0;
  Float_t AngleRcP2 = 0;
  Float_t AngleP1P2 = 0;

  DLM_Random ranGen(11);

  TFile* filePR = TFile::Open(kPRFile, "READ");
  if(!filePR || filePR->IsZombie()){
    std::cout << "ERROR: cannot open " << kPRFile << "\n";
  }
  else{
    TNtuple* tuplePR = dynamic_cast<TNtuple*>(filePR->Get("InfoTuple_ClosePairs"));
    if(!tuplePR){
      std::cout << "ERROR: cannot find InfoTuple_ClosePairs in " << kPRFile << "\n";
    }
    else{
      tuplePR->SetBranchAddress("k_D", &k_D);
      tuplePR->SetBranchAddress("P1", &fP1);
      tuplePR->SetBranchAddress("P2", &fP2);
      tuplePR->SetBranchAddress("M1", &fM1);
      tuplePR->SetBranchAddress("M2", &fM2);
      tuplePR->SetBranchAddress("Tau1", &Tau1);
      tuplePR->SetBranchAddress("Tau2", &Tau2);
      tuplePR->SetBranchAddress("AngleRcP1", &AngleRcP1);
      tuplePR->SetBranchAddress("AngleRcP2", &AngleRcP2);
      tuplePR->SetBranchAddress("AngleP1P2", &AngleP1P2);

      unsigned added = 0;
      for(Long64_t iEntry=0; iEntry<tuplePR->GetEntries(); ++iEntry){
        tuplePR->GetEntry(iEntry);
        if(k_D > cutoffMeV) continue;
        if(fP2 <= 0.f || fM2 <= 0.f || tauReso <= 0.0) continue;

        const double shift = ranGen.Exponential(fM2/(fP2*tauReso));
        source->AddBGT_PR(shift, -std::cos(AngleRcP2));
        source->AddBGT_RP(shift,  std::cos(AngleRcP2));
        ++added;
      }
      std::cout << "Loaded source-shape PR/RP kernels: " << added << "\n";
    }
    filePR->Close();
  }
  delete filePR;

  TFile* fileRR = TFile::Open(kRRFile, "READ");
  if(!fileRR || fileRR->IsZombie()){
    std::cout << "ERROR: cannot open " << kRRFile << "\n";
  }
  else{
    TNtuple* tupleRR = dynamic_cast<TNtuple*>(fileRR->Get("InfoTuple_ClosePairs"));
    if(!tupleRR){
      std::cout << "ERROR: cannot find InfoTuple_ClosePairs in " << kRRFile << "\n";
    }
    else{
      tupleRR->SetBranchAddress("k_D", &k_D);
      tupleRR->SetBranchAddress("P1", &fP1);
      tupleRR->SetBranchAddress("P2", &fP2);
      tupleRR->SetBranchAddress("M1", &fM1);
      tupleRR->SetBranchAddress("M2", &fM2);
      tupleRR->SetBranchAddress("Tau1", &Tau1);
      tupleRR->SetBranchAddress("Tau2", &Tau2);
      tupleRR->SetBranchAddress("AngleRcP1", &AngleRcP1);
      tupleRR->SetBranchAddress("AngleRcP2", &AngleRcP2);
      tupleRR->SetBranchAddress("AngleP1P2", &AngleP1P2);

      unsigned added = 0;
      for(Long64_t iEntry=0; iEntry<tupleRR->GetEntries(); ++iEntry){
        tupleRR->GetEntry(iEntry);
        if(k_D > cutoffMeV) continue;
        if(fM1 <= 0.f || fM2 <= 0.f) continue;
        if(fP1 <= 0.f || fP2 <= 0.f || tauReso <= 0.0) continue;

        const double shift1 = ranGen.Exponential(fM1/(fP1*tauReso));
        const double shift2 = ranGen.Exponential(fM2/(fP2*tauReso));
        source->AddBGT_RR(shift1, std::cos(AngleRcP1),
                          shift2, std::cos(AngleRcP2),
                          std::cos(AngleP1P2));
        ++added;
      }
      std::cout << "Loaded source-shape RR kernels: " << added << "\n";
    }
    fileRR->Close();
  }
  delete fileRR;

  source->InitNumMcIter(1000000);

  double x = 0.1;
  double pars[2] = {radius, alpha};
  source->RootEval(&x, pars);

  return source;
}

void NormalizeGraphArea(TGraph& graph)
{
  double area = 0.0;
  for(int i=0; i<graph.GetN()-1; ++i){
    double x1 = 0.0;
    double y1 = 0.0;
    double x2 = 0.0;
    double y2 = 0.0;
    graph.GetPoint(i, x1, y1);
    graph.GetPoint(i+1, x2, y2);
    area += 0.5*(y1+y2)*(x2-x1);
  }
  if(area <= 0.0) return;
  for(int i=0; i<graph.GetN(); ++i){
    double x = 0.0;
    double y = 0.0;
    graph.GetPoint(i, x, y);
    graph.SetPoint(i, x, y/area);
  }
}

} // namespace

void draw_source_shape_gauss_vs_reso()
{
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  const double rPure = 7.16884;
  const double rReso = 6.88115;
  const double alphaReso = 2.0;
  const double fracReso = 0.72;
  const double tauReso = 1.5;
  const double massReso = 1124.0;
  const double cutoffMeV = 400.0;

  DLM_CleverMcLevyResoTM* reso =
      BuildResonanceSource(rReso, alphaReso, fracReso, tauReso, massReso, cutoffMeV);

  const int n = 600;
  const double rMin = 0.0;
  const double rMax = 80.0;
  TGraph gPureRadial(n);
  TGraph gResoRadial(n);

  double parsReso[2] = {rReso, alphaReso};
  for(int i=0; i<n; ++i){
    const double r = rMin + (rMax-rMin)*double(i)/double(n-1);
    const double sPure = EvalPureGaussSource(r, rPure);
    const double sReso = reso->RootEval((double*)&r, parsReso);

    gPureRadial.SetPoint(i, r, sPure);
    gResoRadial.SetPoint(i, r, sReso);
  }

  NormalizeGraphArea(gPureRadial);
  NormalizeGraphArea(gResoRadial);

  gPureRadial.SetName("gPureGauss_radialPdf");
  gResoRadial.SetName("gGaussResonance_alpha2_radialPdf");
  gPureRadial.SetLineColor(kBlue+1);
  gPureRadial.SetLineWidth(4);
  gResoRadial.SetLineColor(kRed+1);
  gResoRadial.SetLineWidth(4);
  gResoRadial.SetLineStyle(2);

  TCanvas* c1 = new TCanvas("c_source_radial_pdf", "c_source_radial_pdf", 1000, 750);
  gPureRadial.SetTitle(";r* (fm);normalized radial source 4#pi r^{2} S(r*)");
  gPureRadial.GetXaxis()->SetLimits(0, 80);
  gPureRadial.GetYaxis()->SetRangeUser(0, 0.08);
  gPureRadial.Draw("AL");
  gResoRadial.Draw("L SAME");
  TLegend* leg1 = new TLegend(0.48, 0.68, 0.88, 0.86);
  leg1->SetBorderSize(0);
  leg1->SetFillStyle(0);
  leg1->SetTextSize(0.035);
  leg1->AddEntry(&gPureRadial, Form("pure Gauss, R=%.5g fm", rPure), "l");
  leg1->AddEntry(&gResoRadial, Form("Gauss+resonance, R=%.5g fm, #alpha=2", rReso), "l");
  leg1->Draw();

  TCanvas* c2 = new TCanvas("c_source_density", "c_source_density", 1000, 750);
  c2->SetLogy();
  gPureRadial.SetTitle(";r* (fm);normalized radial source 4#pi r^{2} S(r*)");
  gPureRadial.GetXaxis()->SetLimits(0, 80);
  gPureRadial.GetYaxis()->SetRangeUser(1e-6, 0.08);
  gPureRadial.Draw("AL");
  gResoRadial.Draw("L SAME");
  TLegend* leg2 = new TLegend(0.48, 0.68, 0.88, 0.86);
  leg2->SetBorderSize(0);
  leg2->SetFillStyle(0);
  leg2->SetTextSize(0.035);
  leg2->AddEntry(&gPureRadial, Form("pure Gauss, R=%.5g fm", rPure), "l");
  leg2->AddEntry(&gResoRadial, Form("Gauss+resonance, R=%.5g fm, #alpha=2", rReso), "l");
  leg2->Draw();

  const TString outRoot = kOutDir + "/source_shape_pureGauss_vs_resoAlpha2.root";
  const TString outRadialPng = kOutDir + "/source_shape_radialPdf_pureGauss_vs_resoAlpha2.png";
  const TString outDensityPng = kOutDir + "/source_shape_density_pureGauss_vs_resoAlpha2.png";
  const TString outRadialPdf = kOutDir + "/source_shape_radialPdf_pureGauss_vs_resoAlpha2.pdf";
  const TString outDensityPdf = kOutDir + "/source_shape_density_pureGauss_vs_resoAlpha2.pdf";

  c1->SaveAs(outRadialPng);
  c1->SaveAs(outRadialPdf);
  c2->SaveAs(outDensityPng);
  c2->SaveAs(outDensityPdf);

  TFile out(outRoot, "RECREATE");
  gPureRadial.Write();
  gResoRadial.Write();
  c1->Write();
  c2->Write();
  out.Close();

  std::cout << "Saved source-shape plots:\n"
            << "  " << outRadialPng << "\n"
            << "  " << outDensityPng << "\n"
            << "  " << outRoot << "\n";

  delete reso;
}
