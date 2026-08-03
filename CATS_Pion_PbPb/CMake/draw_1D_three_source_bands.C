#include "TCanvas.h"
#include "TFile.h"
#include "TF1.h"
#include "TGraphAsymmErrors.h"
#include "TH1.h"
#include "TLegend.h"
#include "TLine.h"
#include "TPad.h"
#include "TStyle.h"
#include "TString.h"

#include <algorithm>
#include <iostream>
#include <vector>

namespace {

struct FitCurve {
  TString label;
  TString fileName;
  TString dataName;
  TString fitName;
  TString fitFuncName;
  TString bandName;
  Color_t color;
  Color_t fillColor;
  TH1* data = nullptr;
  TH1* fit = nullptr;
  TF1* fitFunc = nullptr;
  TGraphAsymmErrors* band = nullptr;
  TGraphAsymmErrors* ratio = nullptr;
  TGraphAsymmErrors* ratioBand = nullptr;
};

template <class TObjectType>
TObjectType* GetClone(TFile& file, const TString& name, const TString& newName)
{
  TObjectType* obj = dynamic_cast<TObjectType*>(file.Get(name));
  if(!obj){
    std::cout << "ERROR: missing object " << name << " in "
              << file.GetName() << "\n";
    return nullptr;
  }

  TObjectType* copy = dynamic_cast<TObjectType*>(obj->Clone(newName));
  if(auto* hist = dynamic_cast<TH1*>(copy)) hist->SetDirectory(nullptr);
  return copy;
}

TGraphAsymmErrors* BuildDataOverFitRatio(const TH1* data,
                                         const TH1* fit,
                                         const TString& name)
{
  if(!data || !fit) return nullptr;

  TGraphAsymmErrors* ratio = new TGraphAsymmErrors(data->GetNbinsX());
  ratio->SetName(name);
  ratio->SetTitle(name);

  for(int iBin=1; iBin<=data->GetNbinsX(); ++iBin){
    const double x = data->GetXaxis()->GetBinCenter(iBin);
    const double ex = 0.5*data->GetXaxis()->GetBinWidth(iBin);
    const double yData = data->GetBinContent(iBin);
    const double eData = data->GetBinError(iBin);
    const double yFit = fit->GetBinContent(iBin);

    double y = 0.0;
    double ey = 0.0;
    if(yFit > 0.0){
      y = yData/yFit;
      ey = eData/yFit;
    }

    ratio->SetPoint(iBin-1, x, y);
    ratio->SetPointError(iBin-1, ex, ex, ey, ey);
  }

  return ratio;
}

TGraphAsymmErrors* BuildDataOverFitRatio(const TH1* data,
                                         const TF1* fit,
                                         const TString& name)
{
  if(!data || !fit) return nullptr;

  TGraphAsymmErrors* ratio = new TGraphAsymmErrors(data->GetNbinsX());
  ratio->SetName(name);
  ratio->SetTitle(name);

  for(int iBin=1; iBin<=data->GetNbinsX(); ++iBin){
    const double x = data->GetXaxis()->GetBinCenter(iBin);
    const double ex = 0.5*data->GetXaxis()->GetBinWidth(iBin);
    const double yData = data->GetBinContent(iBin);
    const double eData = data->GetBinError(iBin);
    const double yFit = fit->Eval(x);

    double y = 0.0;
    double ey = 0.0;
    if(yFit > 0.0){
      y = yData/yFit;
      ey = eData/yFit;
    }

    ratio->SetPoint(iBin-1, x, y);
    ratio->SetPointError(iBin-1, ex, ex, ey, ey);
  }

  return ratio;
}

TGraphAsymmErrors* BuildRatioBand(const TGraphAsymmErrors* band,
                                  const TH1* fit,
                                  const TString& name)
{
  if(!band || !fit) return nullptr;

  TGraphAsymmErrors* out = new TGraphAsymmErrors(fit->GetNbinsX());
  out->SetName(name);
  out->SetTitle(name);

  for(int iBin=1; iBin<=fit->GetNbinsX(); ++iBin){
    const int i = iBin-1;
    double x = 0.0;
    double y = 0.0;
    band->GetPoint(i, x, y);
    const double yFit = fit->GetBinContent(iBin);
    const double relLow = yFit > 0.0 ? band->GetErrorYlow(i)/yFit : 0.0;
    const double relHigh = yFit > 0.0 ? band->GetErrorYhigh(i)/yFit : 0.0;

    out->SetPoint(i, x, 1.0);
    out->SetPointError(i,
                       band->GetErrorXlow(i),
                       band->GetErrorXhigh(i),
                       relLow,
                       relHigh);
  }

  return out;
}

TGraphAsymmErrors* ClipGraphX(const TGraphAsymmErrors* in,
                              const TString& name,
                              const double xmin,
                              const double xmax)
{
  if(!in) return nullptr;

  TGraphAsymmErrors* out = new TGraphAsymmErrors();
  out->SetName(name);
  out->SetTitle(name);

  int outPoint = 0;
  for(int i=0; i<in->GetN(); ++i){
    double x = 0.0;
    double y = 0.0;
    in->GetPoint(i, x, y);
    if(x < xmin || x > xmax) continue;

    out->SetPoint(outPoint, x, y);
    out->SetPointError(outPoint,
                       in->GetErrorXlow(i),
                       in->GetErrorXhigh(i),
                       in->GetErrorYlow(i),
                       in->GetErrorYhigh(i));
    ++outPoint;
  }

  return out;
}

bool LoadCurve(FitCurve& curve, const TString& tag)
{
  TFile file(curve.fileName, "READ");
  if(file.IsZombie()){
    std::cout << "ERROR: cannot open " << curve.fileName << "\n";
    return false;
  }

  curve.data = GetClone<TH1>(file, curve.dataName, tag + "_data");
  curve.fit = GetClone<TH1>(file, curve.fitName, tag + "_fit");
  curve.fitFunc = GetClone<TF1>(file, curve.fitFuncName, tag + "_fitFunc");
  curve.band = GetClone<TGraphAsymmErrors>(file, curve.bandName, tag + "_band");
  if(!curve.data || !curve.fit || !curve.fitFunc || !curve.band) return false;

  curve.ratio = BuildDataOverFitRatio(curve.data, curve.fitFunc, tag + "_data_over_fit");
  curve.ratioBand = BuildRatioBand(curve.band, curve.fit, tag + "_ratio_band");
  return curve.ratio && curve.ratioBand;
}

void StyleCurve(FitCurve& curve, int markerStyle, int lineStyle, int lineWidth)
{
  curve.fit->SetLineColor(curve.color);
  curve.fit->SetLineWidth(lineWidth);
  curve.fit->SetLineStyle(lineStyle);
  curve.fit->SetMarkerSize(0);
  curve.fitFunc->SetLineColor(curve.color);
  curve.fitFunc->SetLineWidth(lineWidth);
  curve.fitFunc->SetLineStyle(lineStyle);
  curve.fitFunc->SetRange(2.0, 80.0);

  curve.band->SetFillColorAlpha(curve.fillColor, 0.22);
  curve.band->SetLineColorAlpha(curve.fillColor, 0.0);
  curve.band->SetMarkerSize(0);

  curve.ratio->SetMarkerColor(curve.color);
  curve.ratio->SetLineColor(curve.color);
  curve.ratio->SetMarkerStyle(markerStyle);
  curve.ratio->SetMarkerSize(0.85);
  curve.ratio->SetLineWidth(1);

  curve.ratioBand->SetFillColorAlpha(curve.fillColor, 0.16);
  curve.ratioBand->SetLineColorAlpha(curve.fillColor, 0.0);
  curve.ratioBand->SetMarkerSize(0);
}

} // namespace

void draw_1D_three_source_bands()
{
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  const double plotXMin = 2.0;
  const double plotXMax = 80.0;

  const TString baseDir =
      "/Users/wuwenya/Documents/code/femto_flow/CATS/CATS_install/CATS_Pion_PbPb/CMake/out/Aug";

  std::vector<FitCurve> curves = {
      {"pure Gaussian, R=7.169 fm",
       baseDir + "/fit_results_gaussOnly_ForMax_cent10to30_mT4.root",
       "copyhisto_gaussOnly",
       "CATS_fit_gaussOnly_central",
       "CATS_fit_gaussOnly",
       "CATS_fit_gaussOnly_sourceR_1sigma_band",
       kMagenta+1,
       kAzure+1},
      {"Gauss + resonance, #alpha=2.0, R=7.261 fm",
       baseDir + "/fit_results_CleverMcLevyResoTM_ForMax_fixedAlpha2p0_cent10to30_mT4.root",
       "copyhisto_CleverMcLevyResoTM",
       "CATS_fit_CleverMcLevyResoTM_central",
       "CATS_fit_CleverMcLevyResoTM",
       "CATS_fit_CleverMcLevyResoTM_sourceR_1sigma_band",
       kBlue+1,
       kRed-9},
      {"Gauss + resonance, #alpha=1.6, R=8.056 fm",
       baseDir + "/fit_results_CleverMcLevyResoTM_ForMax_fixedAlpha1p6_cent10to30_mT4.root",
       "copyhisto_CleverMcLevyResoTM",
       "CATS_fit_CleverMcLevyResoTM_central",
       "CATS_fit_CleverMcLevyResoTM",
       "CATS_fit_CleverMcLevyResoTM_sourceR_1sigma_band",
       kOrange+7,
       kGreen-7}
  };

  for(size_t i=0; i<curves.size(); ++i){
    if(!LoadCurve(curves[i], Form("curve%zu", i))){
      std::cout << "ERROR: failed to load input curve " << i << "\n";
      return;
    }
    TGraphAsymmErrors* clippedBand =
        ClipGraphX(curves[i].band, Form("curve%zu_band_plot", i), plotXMin, plotXMax);
    TGraphAsymmErrors* clippedRatio =
        ClipGraphX(curves[i].ratio, Form("curve%zu_ratio_plot", i), plotXMin, plotXMax);
    TGraphAsymmErrors* clippedRatioBand =
        ClipGraphX(curves[i].ratioBand, Form("curve%zu_ratio_band_plot", i), plotXMin, plotXMax);
    delete curves[i].band;
    delete curves[i].ratio;
    delete curves[i].ratioBand;
    curves[i].band = clippedBand;
    curves[i].ratio = clippedRatio;
    curves[i].ratioBand = clippedRatioBand;
    const int lineStyle = i==1 ? 2 : (i==2 ? 7 : 1);
    const int lineWidth = i==1 ? 7 : 6;
    StyleCurve(curves[i], 20 + int(i), lineStyle, lineWidth);
  }

  TH1* hData = dynamic_cast<TH1*>(curves.front().data->Clone("hDataForPlot"));
  hData->SetDirectory(nullptr);
  hData->SetMarkerStyle(24);
  hData->SetMarkerSize(1.15);
  hData->SetMarkerColor(kBlack);
  hData->SetLineColor(kBlack);
  hData->SetTitle("");
  hData->GetXaxis()->SetTitle("");
  hData->GetYaxis()->SetTitle("C(#it{k*})");
  hData->GetYaxis()->SetTitleSize(0.058);
  hData->GetYaxis()->SetLabelSize(0.050);
  hData->GetYaxis()->SetTitleOffset(1.05);
  hData->GetXaxis()->SetLabelSize(0.0);
  hData->GetYaxis()->SetRangeUser(0.96, 1.34);
  hData->GetXaxis()->SetRangeUser(plotXMin, plotXMax);

  TCanvas* c = new TCanvas("c_1D_source_compare_bands",
                           "c_1D_source_compare_bands",
                           1200, 900);
  TPad* pTop = new TPad("pTop", "pTop", 0.0, 0.30, 1.0, 1.0);
  TPad* pBot = new TPad("pBot", "pBot", 0.0, 0.0, 1.0, 0.30);
  pTop->SetBottomMargin(0.02);
  pTop->SetLeftMargin(0.15);
  pTop->SetRightMargin(0.04);
  pTop->SetTopMargin(0.04);
  pBot->SetTopMargin(0.02);
  pBot->SetBottomMargin(0.32);
  pBot->SetLeftMargin(0.15);
  pBot->SetRightMargin(0.04);
  pTop->Draw();
  pBot->Draw();

  pTop->cd();
  hData->Draw("E1");
  curves[2].fitFunc->Draw("SAME");
  curves[1].fitFunc->Draw("SAME");
  curves[0].fitFunc->Draw("SAME");
  hData->Draw("E1 SAME");

  TLegend* leg = new TLegend(0.36, 0.45, 0.92, 0.85);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.052);
  leg->AddEntry(hData, "data", "pe");
  for(auto& curve : curves){
    leg->AddEntry(curve.fitFunc, curve.label, "l");
  }
  leg->Draw();

  pBot->cd();
  TH1* frame = dynamic_cast<TH1*>(hData->Clone("ratioFrame"));
  frame->SetDirectory(nullptr);
  frame->Reset("ICES");
  frame->GetYaxis()->SetRangeUser(0.90, 1.10);
  frame->GetXaxis()->SetRangeUser(plotXMin, plotXMax);
  frame->GetYaxis()->SetTitle("data / fit");
  frame->GetXaxis()->SetTitle("#it{k*} (MeV/#it{c})");
  frame->GetYaxis()->SetNdivisions(505);
  frame->GetYaxis()->SetTitleSize(0.105);
  frame->GetYaxis()->SetLabelSize(0.085);
  frame->GetYaxis()->SetTitleOffset(0.62);
  frame->GetXaxis()->SetTitleSize(0.105);
  frame->GetXaxis()->SetLabelSize(0.090);
  frame->GetXaxis()->SetTitleOffset(1.02);
  frame->Draw("AXIS");

  TLine* one = new TLine(plotXMin, 1.0, plotXMax, 1.0);
  one->SetLineColor(kGray+2);
  one->SetLineStyle(2);
  one->Draw("SAME");
  for(auto& curve : curves) curve.ratio->Draw("PZ SAME");

  c->cd();
  const TString outPng = baseDir + "/compare_1D_three_source_lines_ForMax_cent10to30_mT4.png";
  const TString outPdf = baseDir + "/compare_1D_three_source_lines_ForMax_cent10to30_mT4.pdf";
  const TString outRoot = baseDir + "/compare_1D_three_source_lines_ForMax_cent10to30_mT4.root";
  c->SaveAs(outPng);
  c->SaveAs(outPdf);

  TFile out(outRoot, "RECREATE");
  c->Write();
  hData->Write("hData");
  for(auto& curve : curves){
    curve.fit->Write();
    curve.fitFunc->Write();
    curve.ratio->Write();
  }
  out.Close();

  std::cout << "Saved comparison plot:\n"
            << "  " << outPng << "\n"
            << "  " << outPdf << "\n"
            << "  " << outRoot << "\n";
}
