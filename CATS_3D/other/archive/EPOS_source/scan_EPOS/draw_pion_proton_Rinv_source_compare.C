#include <TCanvas.h>
#include <TFile.h>
#include <TF1.h>
#include <TH1.h>
#include <TLegend.h>
#include <TString.h>
#include <TStyle.h>

#include <algorithm>
#include <iostream>

namespace {

TH1D* LoadNormalizedHist(TFile* file, const char* histName, const char* cloneName)
{
  TH1D* inputHist = dynamic_cast<TH1D*>(file->Get(histName));
  if (!inputHist) {
    std::cerr << "ERROR: missing histogram " << histName << std::endl;
    return nullptr;
  }

  TH1D* hist = dynamic_cast<TH1D*>(inputHist->Clone(cloneName));
  hist->SetDirectory(nullptr);
  const double integral = hist->Integral("width");
  if (integral > 0.) hist->Scale(1. / integral);
  return hist;
}

TF1* LoadFit(TFile* file, const char* functionName, int color, int lineStyle, int lineWidth)
{
  TF1* inputFunction = dynamic_cast<TF1*>(file->Get(functionName));
  if (!inputFunction) {
    std::cerr << "ERROR: missing fit function " << functionName << std::endl;
    return nullptr;
  }

  TF1* function = dynamic_cast<TF1*>(inputFunction->Clone(Form("%s_draw", functionName)));
  function->SetLineColor(color);
  function->SetLineStyle(lineStyle);
  function->SetLineWidth(lineWidth);
  function->SetNpx(1000);
  return function;
}

}  // namespace

void draw_pion_proton_Rinv_source_compare(
    TString pionFileName = "epos_primordial_pion_pair_Rinv_fit.root",
    TString protonFileName = "epos_primordial_protn_pair_Rinv_fit.root",
    TString outputPrefix = "epos_primordial_pion_proton_Rinv_source_compare")
{
  TFile* pionFile = TFile::Open(pionFileName, "READ");
  TFile* protonFile = TFile::Open(protonFileName, "READ");
  if (!pionFile || pionFile->IsZombie()) {
    std::cerr << "ERROR: cannot open " << pionFileName << std::endl;
    return;
  }
  if (!protonFile || protonFile->IsZombie()) {
    std::cerr << "ERROR: cannot open " << protonFileName << std::endl;
    return;
  }

  TH1D* hPion = LoadNormalizedHist(pionFile, "hRInv_PRF_levy_fit_input", "hPionRInvNorm");
  TH1D* hProton = LoadNormalizedHist(protonFile, "hRInv_PRF_levy_fit_input", "hProtonRInvNorm");
  TF1* fPionLogNormal = LoadFit(pionFile, "normLogNormalRadiusMixture", kRed + 1, 1, 3);
  TF1* fPionQGaussian = LoadFit(pionFile, "normGammaMixtureQGaussian", kMagenta + 2, 2, 3);
  TF1* fProtonGaussian = LoadFit(protonFile, "normPureGaussian", kGreen + 2, 1, 3);

  pionFile->Close();
  protonFile->Close();

  if (!hPion || !hProton || !fPionLogNormal || !fPionQGaussian || !fProtonGaussian) return;

  hPion->SetMarkerStyle(20);
  hPion->SetMarkerSize(0.75);
  hPion->SetMarkerColor(kBlack);
  hPion->SetLineColor(kBlack);
  hPion->SetLineWidth(1);

  hProton->SetMarkerStyle(24);
  hProton->SetMarkerSize(0.85);
  hProton->SetMarkerColor(kBlue + 1);
  hProton->SetLineColor(kBlue + 1);
  hProton->SetLineWidth(1);

  const double xMin = 0.;
  const double xMax = 200.;
  hPion->GetXaxis()->SetRangeUser(xMin, xMax);
  hProton->GetXaxis()->SetRangeUser(xMin, xMax);

  double yMax = 0.;
  const int pionFirst = hPion->GetXaxis()->FindBin(xMin);
  const int pionLast = hPion->GetXaxis()->FindBin(xMax);
  for (int iBin = pionFirst; iBin <= pionLast; ++iBin) {
    yMax = std::max(yMax, hPion->GetBinContent(iBin) + hPion->GetBinError(iBin));
  }
  const int protonFirst = hProton->GetXaxis()->FindBin(xMin);
  const int protonLast = hProton->GetXaxis()->FindBin(xMax);
  for (int iBin = protonFirst; iBin <= protonLast; ++iBin) {
    yMax = std::max(yMax, hProton->GetBinContent(iBin) + hProton->GetBinError(iBin));
  }

  gStyle->SetOptStat(0);
  TCanvas* canvas = new TCanvas("cPionProtonRInvSourceCompare",
                                "EPOS pion/proton primordial pair source comparison",
                                950, 720);
  canvas->SetTicks(1, 1);
  canvas->SetGrid(1, 1);

  hPion->SetTitle(";#it{r}_{inv} (fm);4#pi#it{r}^{2}S(#it{r}) (Norm.)");
  hPion->GetYaxis()->SetRangeUser(0., 1.25 * yMax);
  hPion->GetXaxis()->SetTitleOffset(1.05);
  hPion->GetYaxis()->SetTitleOffset(1.25);
  hPion->Draw("E1");
  hProton->Draw("E1 SAME");
  fPionLogNormal->Draw("SAME");
  fPionQGaussian->Draw("SAME");
  fProtonGaussian->Draw("SAME");

  TLegend* legend = new TLegend(0.53, 0.60, 0.88, 0.88);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->SetTextSize(0.032);
  legend->AddEntry(hPion, "EPOS #pi#pi source", "lep");
  legend->AddEntry(fPionLogNormal, "#pi#pi log-normal mixture", "l");
  legend->AddEntry(fPionQGaussian, "#pi#pi q-Gaussian mixture", "l");
  legend->AddEntry(hProton, "EPOS pp source", "lep");
  legend->AddEntry(fProtonGaussian, "pp pure Gaussian", "l");
  legend->Draw();

  canvas->SaveAs(outputPrefix + ".pdf");
  canvas->SaveAs(outputPrefix + ".png");

  TFile outputFile(outputPrefix + ".root", "RECREATE");
  canvas->Write();
  hPion->Write();
  hProton->Write();
  fPionLogNormal->Write();
  fPionQGaussian->Write();
  fProtonGaussian->Write();
  outputFile.Close();

  std::cout << "Saved " << outputPrefix << ".pdf/.png/.root" << std::endl;
}
