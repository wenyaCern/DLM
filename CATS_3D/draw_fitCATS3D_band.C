#include "TCanvas.h"
#include "TFile.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TLine.h"
#include "TString.h"
#include "TStyle.h"
#include "TPad.h"
#include "TAxis.h"

#include <cmath>
#include <iostream>
#include <memory>

namespace {

TH1D* GetHist1D(TFile& file, const TString& name)
{
    TH1D* h = dynamic_cast<TH1D*>(file.Get(name));
    if(!h){
        std::cout << "ERROR: missing histogram " << name << "\n";
        return nullptr;
    }
    h->SetDirectory(nullptr);
    return h;
}

TH1D* GetOptionalHist1D(TFile& file, const TString& name)
{
    TH1D* h = dynamic_cast<TH1D*>(file.Get(name));
    if(!h) return nullptr;
    h->SetDirectory(nullptr);
    return h;
}

TH1D* MakeBandHist(const TH1D* hLow,
                   const TH1D* hHigh,
                   const TString& name)
{
    if(!hLow || !hHigh) return nullptr;

    TH1D* hBand = dynamic_cast<TH1D*>(hLow->Clone(name));
    hBand->SetDirectory(nullptr);

    for(int i=1; i<=hLow->GetNbinsX(); ++i){
        const double yLow = hLow->GetBinContent(i);
        const double yHigh = hHigh->GetBinContent(i);
        const double y = 0.5 * (yLow + yHigh);
        hBand->SetBinContent(i, y);
        hBand->SetBinError(i, 0.5 * std::fabs(yHigh - yLow));
    }

    hBand->SetFillColorAlpha(kAzure-9, 0.65);
    hBand->SetLineColor(kAzure-9);
    hBand->SetMarkerSize(0);
    return hBand;
}

TH1D* MakeRatioHist(const TH1D* hData,
                    const TH1D* hFit,
                    const TString& name)
{
    if(!hData || !hFit) return nullptr;

    TH1D* hRatio = dynamic_cast<TH1D*>(hData->Clone(name));
    hRatio->SetDirectory(nullptr);
    hRatio->Reset("ICES");

    for(int i=1; i<=hData->GetNbinsX(); ++i){
        const double d = hData->GetBinContent(i);
        const double ed = hData->GetBinError(i);
        const double f = hFit->GetBinContent(i);
        if(f == 0.0) continue;
        hRatio->SetBinContent(i, d/f);
        hRatio->SetBinError(i, ed/std::fabs(f));
    }
    return hRatio;
}

TH1D* MakeRatioBandHist(const TH1D* hLow,
                        const TH1D* hHigh,
                        const TH1D* hFit,
                        const TString& name)
{
    if(!hLow || !hHigh || !hFit) return nullptr;

    TH1D* hBand = dynamic_cast<TH1D*>(hFit->Clone(name));
    hBand->SetDirectory(nullptr);
    hBand->Reset("ICES");

    for(int i=1; i<=hLow->GetNbinsX(); ++i){
        const double f = hFit->GetBinContent(i);
        if(f == 0.0){
            hBand->SetBinContent(i, 1.0);
            hBand->SetBinError(i, 0.0);
            continue;
        }
        const double yLow = hLow->GetBinContent(i) / f;
        const double yHigh = hHigh->GetBinContent(i) / f;
        const double y = 0.5 * (yLow + yHigh);
        hBand->SetBinContent(i, y);
        hBand->SetBinError(i, 0.5 * std::fabs(yHigh - yLow));
    }

    hBand->SetFillColorAlpha(kAzure-9, 0.65);
    hBand->SetLineColor(kAzure-9);
    hBand->SetMarkerSize(0);
    return hBand;
}

TString AxisTitle(const TString& direction)
{
    if(direction == "out") return "q_{out} (GeV/c)";
    if(direction == "side") return "q_{side} (GeV/c)";
    return "q_{long} (GeV/c)";
}

} // namespace

void draw_fitCATS3D_band(const char* fileName = "../fit_CATS3D_resonance_cent1_mt4_wtCoulomb_levyAlpha2.0.root",
                         const char* direction = "out",
                         bool useMixedWeighted = true,
                         bool drawBackground = true,
                         const char* outName = "")
{
    gStyle->SetOptStat(0);
    gStyle->SetEndErrorSize(4);

    TFile file(fileName, "READ");
    if(file.IsZombie()){
        std::cout << "ERROR: cannot open " << fileName << "\n";
        return;
    }

    const TString dir(direction);
    const TString projTag = useMixedWeighted ? "_mixedWeighted_proj_" : "_proj_";
    const TString bandTag = useMixedWeighted
        ? "hCATS3D_best_1sigma_mixedWeighted_proj_"
        : "hCATS3D_best_1sigma_proj_";
    TH1D* hData = GetHist1D(file, "hDataCk3D" + projTag + dir);
    TH1D* hFit = GetHist1D(file, "hCATS3D_best" + projTag + dir);
    TH1D* hLow = GetHist1D(file, bandTag + dir + "_low");
    TH1D* hHigh = GetHist1D(file, bandTag + dir + "_high");
    if(!hData || !hFit || !hLow || !hHigh) return;

    TH1D* hBkg = nullptr;
    TH1D* hBkgBand = nullptr;
    if(drawBackground){
        hBkg = GetOptionalHist1D(file, "hCATS3D_background" + projTag + dir);
        TH1D* hBkgLow = GetOptionalHist1D(file,
            "hCATS3D_background_1sigma_low" + projTag + dir);
        TH1D* hBkgHigh = GetOptionalHist1D(file,
            "hCATS3D_background_1sigma_high" + projTag + dir);
        if(hBkgLow && hBkgHigh){
            hBkgBand = MakeBandHist(hBkgLow, hBkgHigh,
                                    "hBkgBand_proj_" + dir);
            hBkgBand->SetFillColorAlpha(kOrange+1, 0.45);
            hBkgBand->SetLineColor(kOrange+7);
            hBkgBand->SetLineWidth(1);
        }
        if(!hBkg){
            std::cout << "WARNING: missing background projection "
                      << "hCATS3D_background" + projTag + dir << "\n";
        }
        if(!hBkgBand){
            std::cout << "WARNING: missing background 1-sigma band. "
                      << "Rerun the fit after the bkg-band code update.\n";
        }
    }

    TH1D* hRatio = MakeRatioHist(hData, hFit, "hDataOverFit_proj_" + dir);
    TH1D* hBand = MakeBandHist(
        hLow, hHigh, "hFitBand_proj_" + dir);
    TH1D* hRatioBand = MakeRatioBandHist(
        hLow, hHigh, hFit, "hFitRatioBand_proj_" + dir);

    hData->SetMarkerStyle(20);
    hData->SetMarkerSize(1.15);
    hData->SetMarkerColor(kBlack);
    hData->SetLineColor(kBlack);

    hFit->SetLineColor(kBlue+1);
    hFit->SetLineWidth(2);
    hFit->SetMarkerSize(0);

    if(hBkg){
        hBkg->SetLineColor(kOrange+7);
        hBkg->SetLineWidth(2);
        hBkg->SetLineStyle(2);
        hBkg->SetMarkerSize(0);
    }

    hRatio->SetMarkerStyle(20);
    hRatio->SetMarkerSize(1.0);
    hRatio->SetMarkerColor(kBlack);
    hRatio->SetLineColor(kBlack);

    auto* c = new TCanvas(Form("cFitBand_%s", dir.Data()),
                          Form("CATS3D fit band %s", dir.Data()),
                          900, 850);

    auto* padTop = new TPad("padTop", "padTop", 0.0, 0.30, 1.0, 1.0);
    auto* padBot = new TPad("padBot", "padBot", 0.0, 0.0, 1.0, 0.30);
    padTop->SetBottomMargin(0.02);
    padTop->SetLeftMargin(0.13);
    padTop->SetRightMargin(0.04);
    padBot->SetTopMargin(0.03);
    padBot->SetBottomMargin(0.32);
    padBot->SetLeftMargin(0.13);
    padBot->SetRightMargin(0.04);
    padTop->Draw();
    padBot->Draw();

    padTop->cd();
    hData->SetTitle("");
    hData->GetXaxis()->SetLabelSize(0);
    hData->GetXaxis()->SetTitleSize(0);
    hData->GetYaxis()->SetTitle("C(q)");
    hData->GetYaxis()->SetTitleSize(0.065);
    hData->GetYaxis()->SetTitleOffset(0.75);
    hData->GetYaxis()->SetLabelSize(0.052);

    double yMin = hData->GetMinimum();
    double yMax = hData->GetMaximum();
    for(int i=1; i<=hLow->GetNbinsX(); ++i){
        yMin = std::min(yMin, hLow->GetBinContent(i));
        yMax = std::max(yMax, hHigh->GetBinContent(i));
    }
    if(hBkgBand){
        for(int i=1; i<=hBkgBand->GetNbinsX(); ++i){
            yMin = std::min(yMin, hBkgBand->GetBinContent(i) -
                                  hBkgBand->GetBinError(i));
            yMax = std::max(yMax, hBkgBand->GetBinContent(i) +
                                  hBkgBand->GetBinError(i));
        }
    }
    hData->GetYaxis()->SetRangeUser(yMin - 0.08*(yMax-yMin),
                                    yMax + 0.22*(yMax-yMin));
    hData->Draw("E1");
    hBand->Draw("E2 SAME");
    if(hBkgBand) hBkgBand->Draw("E2 SAME");
    if(hBkg) hBkg->Draw("HIST SAME");
    hFit->Draw("HIST SAME");
    hData->Draw("E1 SAME");

    auto* leg = new TLegend(0.62, 0.66, 0.93, 0.90);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.048);
    leg->AddEntry(hData, "data", "lep");
    leg->AddEntry(hFit, "CATS3D fit", "l");
    leg->AddEntry(hBand, "1#sigma source band", "f");
    if(hBkg) leg->AddEntry(hBkg, "background", "l");
    if(hBkgBand) leg->AddEntry(hBkgBand, "1#sigma bkg band", "f");
    leg->Draw();

    padBot->cd();
    hRatio->SetTitle("");
    hRatio->GetXaxis()->SetTitle(AxisTitle(dir));
    hRatio->GetXaxis()->SetTitleSize(0.12);
    hRatio->GetXaxis()->SetTitleOffset(1.0);
    hRatio->GetXaxis()->SetLabelSize(0.10);
    hRatio->GetYaxis()->SetTitle("data / fit");
    hRatio->GetYaxis()->SetTitleSize(0.10);
    hRatio->GetYaxis()->SetTitleOffset(0.42);
    hRatio->GetYaxis()->SetLabelSize(0.085);
    hRatio->GetYaxis()->SetNdivisions(505);
    hRatio->GetYaxis()->SetRangeUser(0.92, 1.08);
    hRatio->Draw("E1");
    hRatioBand->Draw("E2 SAME");
    hRatio->Draw("E1 SAME");

    auto* line = new TLine(hRatio->GetXaxis()->GetXmin(), 1.0,
                           hRatio->GetXaxis()->GetXmax(), 1.0);
    line->SetLineColor(kGray+2);
    line->SetLineStyle(2);
    line->Draw("SAME");

    c->cd();
    c->Update();

    TString saveName(outName);
    if(saveName.Length() == 0){
        saveName = Form("fitCATS3D_band_%s.png", dir.Data());
    }
    c->SaveAs(saveName);
}
