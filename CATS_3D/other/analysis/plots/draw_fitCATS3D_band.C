#include "TCanvas.h"
#include "TFile.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TLine.h"
#include "TString.h"
#include "TStyle.h"
#include "TPad.h"
#include "TAxis.h"
#include "TNtuple.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <memory>
#include <vector>

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

double SubSampleP0Error(const char* pattern,
                        const int nSubSamples,
                        const bool divideBySqrtN)
{
    std::vector<double> values;
    for(int i=0; i<nSubSamples; ++i){
        TFile f(Form(pattern, i), "READ");
        TNtuple* t = f.IsZombie() ? nullptr
                     : dynamic_cast<TNtuple*>(f.Get("fit_summary"));
        if(!t || !t->GetLeaf("BkgP0") || t->GetEntries()<1) continue;
        t->GetEntry(0);
        values.push_back(t->GetLeaf("BkgP0")->GetValue());
    }
    if(values.size()<2) return 0.0;
    double mean=0.0;
    for(double v: values) mean += v;
    mean /= values.size();
    double variance=0.0;
    for(double v: values) variance += (v-mean)*(v-mean);
    variance /= double(values.size()-1);
    double error=std::sqrt(variance);
    if(divideBySqrtN) error /= std::sqrt(double(values.size()));
    std::cout << "baseline p0 statistical error=" << error
              << " from " << values.size() << " subSamples"
              << (divideBySqrtN ? " (RMS/sqrt(N))" : " (bootstrap RMS)")
              << "\n";
    return error;
}

bool BuildSubSampleFitEnvelope(const TH1D* defaultFit,
                               const char* pattern,
                               const TString& histName,
                               const int nSubSamples,
                               const bool divideBySqrtN,
                               TH1D*& low,
                               TH1D*& high)
{
    if(!defaultFit) return false;
    const int nBins = defaultFit->GetNbinsX();
    std::vector<double> sum(nBins+1, 0.0);
    std::vector<double> sum2(nBins+1, 0.0);
    int nUsed = 0;
    for(int i=0; i<nSubSamples; ++i){
        TFile f(Form(pattern, i), "READ");
        TH1D* h = f.IsZombie() ? nullptr
                  : dynamic_cast<TH1D*>(f.Get(histName));
        if(!h || h->GetNbinsX()!=nBins) continue;
        ++nUsed;
        for(int iBin=1; iBin<=nBins; ++iBin){
            const double value = h->GetBinContent(iBin);
            sum[iBin] += value;
            sum2[iBin] += value*value;
        }
    }
    if(nUsed<2) return false;

    low = dynamic_cast<TH1D*>(defaultFit->Clone("hSubSampleFitLow"));
    high = dynamic_cast<TH1D*>(defaultFit->Clone("hSubSampleFitHigh"));
    low->SetDirectory(nullptr);
    high->SetDirectory(nullptr);
    for(int iBin=1; iBin<=nBins; ++iBin){
        const double mean = sum[iBin]/double(nUsed);
        const double variance = std::max(0.0,
            (sum2[iBin]-double(nUsed)*mean*mean)/double(nUsed-1));
        double error = std::sqrt(variance);
        if(divideBySqrtN) error /= std::sqrt(double(nUsed));
        low->SetBinContent(iBin, defaultFit->GetBinContent(iBin)-error);
        high->SetBinContent(iBin, defaultFit->GetBinContent(iBin)+error);
    }
    std::cout << histName << " statistical band from " << nUsed
              << " bootstrap projections"
              << (divideBySqrtN ? " (bin RMS/sqrt(N))" : " (bootstrap bin RMS)")
              << "\n";
    return true;
}

TH1D* BuildBootstrapResidualSigma(const TH1D* reference,
                                  const char* pattern,
                                  const TString& dataName,
                                  const TString& fitName,
                                  const int nSubSamples,
                                  const bool divideBySqrtN)
{
    if(!reference) return nullptr;
    const int nBins = reference->GetNbinsX();
    std::vector<double> sum(nBins+1, 0.0);
    std::vector<double> sum2(nBins+1, 0.0);
    int nUsed = 0;
    for(int i=0; i<nSubSamples; ++i){
        TFile f(Form(pattern, i), "READ");
        auto* data = f.IsZombie() ? nullptr
                     : dynamic_cast<TH1D*>(f.Get(dataName));
        auto* fit = f.IsZombie() ? nullptr
                    : dynamic_cast<TH1D*>(f.Get(fitName));
        if(!data || !fit || data->GetNbinsX()!=nBins ||
           fit->GetNbinsX()!=nBins) continue;
        ++nUsed;
        for(int iBin=1; iBin<=nBins; ++iBin){
            const double residual = data->GetBinContent(iBin)
                                  - fit->GetBinContent(iBin);
            sum[iBin] += residual;
            sum2[iBin] += residual*residual;
        }
    }
    if(nUsed<2) return nullptr;

    auto* sigmaHist = dynamic_cast<TH1D*>(reference->Clone(
        "hBootstrapResidualSigma"));
    sigmaHist->SetDirectory(nullptr);
    sigmaHist->Reset("ICES");
    for(int iBin=1; iBin<=nBins; ++iBin){
        const double mean = sum[iBin]/double(nUsed);
        const double variance = std::max(0.0,
            (sum2[iBin]-double(nUsed)*mean*mean)/double(nUsed-1));
        double sigma = std::sqrt(variance);
        if(divideBySqrtN) sigma /= std::sqrt(double(nUsed));
        sigmaHist->SetBinError(iBin, sigma);
    }
    std::cout << "paired data-fit residual sigma from " << nUsed
              << " bootstrap projections"
              << (divideBySqrtN ? " (RMS/sqrt(N))" : " (bootstrap RMS)")
              << "\n";
    return sigmaHist;
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

TH1D* MakePullHist(const TH1D* hData,
                   const TH1D* hFit,
                   const TH1D* hBootstrapBand,
                   const TString& name)
{
    if(!hData || !hFit) return nullptr;
    TH1D* hPull = dynamic_cast<TH1D*>(hData->Clone(name));
    hPull->SetDirectory(nullptr);
    hPull->Reset("ICES");
    for(int i=1; i<=hData->GetNbinsX(); ++i){
        const double sigma = hBootstrapBand && hBootstrapBand->GetBinError(i)>0.0
            ? hBootstrapBand->GetBinError(i)
            : hData->GetBinError(i);
        if(!(sigma>0.0) || !std::isfinite(sigma)) continue;
        hPull->SetBinContent(i,
            (hData->GetBinContent(i)-hFit->GetBinContent(i))/sigma);
        hPull->SetBinError(i, 0.0);
    }
    return hPull;
}

TString AxisTitle(const TString& direction)
{
    if(direction == "out") return "q_{out} (GeV/c)";
    if(direction == "side") return "q_{side} (GeV/c)";
    return "q_{long} (GeV/c)";
}

} // namespace

void draw_fitCATS3D_band(const char* fileName = "./fit_CATS3D_gaussian_cent1_mt7_wtCoulomb.root",
                         const char* direction = "side",
                         bool useMixedWeighted = false,
                         bool drawBackground = true,
                         const char* outName = "",
                         const char* subSamplePattern = "",
                         int nSubSamples = 100,
                         bool divideBySqrtN = false,
                         double plotYMin = 0.0,
                         double plotYMax = 0.0)
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
    TH1D* hLow = GetOptionalHist1D(file, bandTag + dir + "_low");
    TH1D* hHigh = GetOptionalHist1D(file, bandTag + dir + "_high");
    if(!hData || !hFit) return;
    TH1D* hDataLow = nullptr;
    TH1D* hDataHigh = nullptr;
    if(subSamplePattern && TString(subSamplePattern).Length()>0){
        BuildSubSampleFitEnvelope(hData, subSamplePattern,
                                  "hDataCk3D" + projTag + dir,
                                  nSubSamples, divideBySqrtN,
                                  hDataLow, hDataHigh);
    }
    if(subSamplePattern && TString(subSamplePattern).Length()>0){
        TH1D* subLow = nullptr;
        TH1D* subHigh = nullptr;
        // Use the CATS shape (before multiplying by p0/background), so this
        // band propagates the fitted Rout/Rside/Rlong fluctuations only.
        // The p0 uncertainty is shown separately by the baseline band.
        const TString subFitName = "hCATS3D_shape" + projTag + dir;
        if(BuildSubSampleFitEnvelope(hFit, subSamplePattern, subFitName,
                                     nSubSamples, divideBySqrtN,
                                     subLow, subHigh)){
            delete hLow;
            delete hHigh;
            hLow = subLow;
            hHigh = subHigh;
        }
    }
    if(!hLow || !hHigh){
        std::cout << "WARNING: missing 1-sigma fit band; "
                  << "draw data and best fit without the band.\n";
        hLow = nullptr;
        hHigh = nullptr;
    }

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
        if(subSamplePattern && TString(subSamplePattern).Length()>0 && hBkg){
            const double p0Error = SubSampleP0Error(subSamplePattern,
                                                    nSubSamples,
                                                    divideBySqrtN);
            if(p0Error>0.0){
                std::unique_ptr<TH1D> low(dynamic_cast<TH1D*>(hBkg->Clone("hBkgP0Low")));
                std::unique_ptr<TH1D> high(dynamic_cast<TH1D*>(hBkg->Clone("hBkgP0High")));
                for(int i=1; i<=hBkg->GetNbinsX(); ++i){
                    low->SetBinContent(i, hBkg->GetBinContent(i)-p0Error);
                    high->SetBinContent(i, hBkg->GetBinContent(i)+p0Error);
                }
                delete hBkgBand;
                hBkgBand = MakeBandHist(low.get(), high.get(),
                                        "hBkgP0StatBand_proj_"+dir);
                hBkgBand->SetFillColorAlpha(kOrange+1, 0.45);
                hBkgBand->SetLineColor(kOrange+7);
                hBkgBand->SetLineWidth(1);
            }
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
    TH1D* hBand = (hLow && hHigh)
        ? MakeBandHist(hLow, hHigh, "hFitBand_proj_" + dir)
        : nullptr;
    TH1D* hRatioBand = (hLow && hHigh)
        ? MakeRatioBandHist(hLow, hHigh, hFit,
                            "hFitRatioBand_proj_" + dir)
        : nullptr;
    TH1D* hDataBand = (hDataLow && hDataHigh)
        ? MakeBandHist(hDataLow, hDataHigh, "hDataBootstrapBand_proj_" + dir)
        : nullptr;
    TH1D* hDataRatioBand = (hDataLow && hDataHigh)
        ? MakeRatioBandHist(hDataLow, hDataHigh, hFit,
                            "hDataBootstrapRatioBand_proj_" + dir)
        : nullptr;
    if(hDataBand){
        hDataBand->SetFillColorAlpha(kGray+1, 0.35);
        hDataBand->SetLineColor(kGray+1);
        // Present the bootstrap q-bin uncertainty as statistical error bars
        // on the full-data points.  A filled band can be mistaken for a
        // systematic uncertainty.
        for(int i=1; i<=hData->GetNbinsX(); ++i){
            hData->SetBinError(i, hDataBand->GetBinError(i));
        }
    }
    if(hDataRatioBand){
        hDataRatioBand->SetFillColorAlpha(kGray+1, 0.35);
        hDataRatioBand->SetLineColor(kGray+1);
    }

    std::unique_ptr<TH1D> hResidualSigma;
    if(subSamplePattern && TString(subSamplePattern).Length()>0){
        hResidualSigma.reset(BuildBootstrapResidualSigma(
            hData, subSamplePattern,
            "hDataCk3D" + projTag + dir,
            "hCATS3D_best" + projTag + dir,
            nSubSamples, divideBySqrtN));
    }

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

    std::unique_ptr<TH1D> hPull(MakePullHist(
        hData, hFit, hResidualSigma ? hResidualSigma.get() : hDataBand,
        "hDataFitPull_proj_" + dir));
    hPull->SetMarkerStyle(20);
    hPull->SetMarkerSize(1.0);
    hPull->SetMarkerColor(kBlack);
    hPull->SetLineColor(kBlack);

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
    if(hLow && hHigh){
        for(int i=1; i<=hLow->GetNbinsX(); ++i){
            yMin = std::min(yMin, hLow->GetBinContent(i));
            yMax = std::max(yMax, hHigh->GetBinContent(i));
        }
    }
    if(hBkgBand){
        for(int i=1; i<=hBkgBand->GetNbinsX(); ++i){
            yMin = std::min(yMin, hBkgBand->GetBinContent(i) -
                                  hBkgBand->GetBinError(i));
            yMax = std::max(yMax, hBkgBand->GetBinContent(i) +
                                  hBkgBand->GetBinError(i));
        }
    }
    if(plotYMax > plotYMin){
        hData->GetYaxis()->SetRangeUser(plotYMin, plotYMax);
    }
    else {
        hData->GetYaxis()->SetRangeUser(yMin - 0.08*(yMax-yMin),
                                        yMax + 0.22*(yMax-yMin));
    }
    hData->Draw("E1");
    if(hBand) hBand->Draw("E2 SAME");
    if(hBkgBand) hBkgBand->Draw("E2 SAME");
    if(hBkg) hBkg->Draw("HIST SAME");
    hFit->Draw("HIST SAME");
    hData->Draw("E1 SAME");

    auto* leg = new TLegend(0.60, 0.58, 0.93, 0.90);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.048);
    leg->AddEntry(hData, hDataBand ? "data (bootstrap stat.)" : "data", "lep");
    leg->AddEntry(hFit, "CATS3D fit", "l");
    if(hBand) leg->AddEntry(hBand, "1#sigma source band", "f");
    if(hBkg) leg->AddEntry(hBkg, "background", "l");
    if(hBkgBand) leg->AddEntry(hBkgBand, "1#sigma bkg band", "f");
    leg->Draw();

    padBot->cd();
    hPull->SetTitle("");
    hPull->GetXaxis()->SetTitle(AxisTitle(dir));
    hPull->GetXaxis()->SetTitleSize(0.12);
    hPull->GetXaxis()->SetTitleOffset(1.0);
    hPull->GetXaxis()->SetLabelSize(0.10);
    hPull->GetYaxis()->SetTitle("(data-fit)/#sigma_{#Delta,boot}");
    hPull->GetYaxis()->SetTitleSize(0.085);
    hPull->GetYaxis()->SetTitleOffset(0.50);
    hPull->GetYaxis()->SetLabelSize(0.085);
    hPull->GetYaxis()->SetNdivisions(505);
    double maxAbsPull = 3.0;
    for(int i=1; i<=hPull->GetNbinsX(); ++i){
        maxAbsPull = std::max(maxAbsPull, std::fabs(hPull->GetBinContent(i)));
    }
    maxAbsPull = std::min(12.0, 1.15*maxAbsPull);
    hPull->GetYaxis()->SetRangeUser(-maxAbsPull, maxAbsPull);
    hPull->Draw("P");

    auto* line = new TLine(hPull->GetXaxis()->GetXmin(), 0.0,
                           hPull->GetXaxis()->GetXmax(), 0.0);
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
