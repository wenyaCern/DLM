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

namespace {

struct FitSummary {
    double Rout = 0.0;
    double Rside = 0.0;
    double Rlong = 0.0;
    double Alpha = 2.0;
};

TH1D* GetHist1D(TFile& file, const TString& name)
{
    TH1D* h = dynamic_cast<TH1D*>(file.Get(name));
    if(!h){
        std::cout << "ERROR: missing histogram " << name
                  << " in " << file.GetName() << "\n";
        return nullptr;
    }
    TH1D* copy = dynamic_cast<TH1D*>(h->Clone(name + "_copy"));
    copy->SetDirectory(nullptr);
    return copy;
}

FitSummary ReadFitSummary(TFile& file)
{
    FitSummary s;
    TNtuple* t = dynamic_cast<TNtuple*>(file.Get("fit_summary"));
    if(!t){
        std::cout << "WARNING: missing fit_summary in "
                  << file.GetName() << "\n";
        return s;
    }

    Float_t Rout = 0.0f;
    Float_t Rside = 0.0f;
    Float_t Rlong = 0.0f;
    Float_t Alpha = 2.0f;
    t->SetBranchAddress("Rout", &Rout);
    t->SetBranchAddress("Rside", &Rside);
    t->SetBranchAddress("Rlong", &Rlong);
    t->SetBranchAddress("Alpha", &Alpha);
    t->GetEntry(0);

    s.Rout = Rout;
    s.Rside = Rside;
    s.Rlong = Rlong;
    s.Alpha = Alpha;
    return s;
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

TString AxisTitle(const TString& direction)
{
    if(direction == "out") return "q_{out} (GeV/c)";
    if(direction == "side") return "q_{side} (GeV/c)";
    return "q_{long} (GeV/c)";
}

TString DirectionTitle(const TString& direction)
{
    if(direction == "out") return "q_{out}";
    if(direction == "side") return "q_{side}";
    return "q_{long}";
}

TString RadiusLabel(const FitSummary& s)
{
    return Form("R=(%.3f, %.3f, %.3f) fm", s.Rout, s.Rside, s.Rlong);
}

} // namespace

void draw_compare_CATS3D_gauss_vs_reso(
    const char* direction = "out",
    bool useMixedWeighted = true,
    const char* gaussianFileName =
        "fit_CATS3D_gaussian_cent1_mt4_wtCoulomb.root",
    const char* resonanceFileName =
        "fit_CATS3D_resonance_cent1_mt4_wtCoulomb_levyAlpha2.0.root",
    const char* outName = "")
{
    gStyle->SetOptStat(0);
    gStyle->SetEndErrorSize(4);

    TFile fGauss(gaussianFileName, "READ");
    TFile fReso(resonanceFileName, "READ");
    if(fGauss.IsZombie() || fReso.IsZombie()){
        std::cout << "ERROR: cannot open input fit root files.\n";
        return;
    }

    const TString dir(direction);
    const TString projTag = useMixedWeighted ? "_mixedWeighted_proj_" : "_proj_";

    TH1D* hData = GetHist1D(fReso, "hDataCk3D" + projTag + dir);
    TH1D* hGauss = GetHist1D(fGauss, "hCATS3D_best" + projTag + dir);
    TH1D* hReso = GetHist1D(fReso, "hCATS3D_best" + projTag + dir);
    if(!hData || !hGauss || !hReso) return;

    FitSummary gaussSummary = ReadFitSummary(fGauss);
    FitSummary resoSummary = ReadFitSummary(fReso);

    TH1D* hRatioGauss = MakeRatioHist(hData, hGauss,
                                      "hDataOverGaussFit_" + dir);
    TH1D* hRatioReso = MakeRatioHist(hData, hReso,
                                     "hDataOverResoFit_" + dir);
    if(!hRatioGauss || !hRatioReso) return;

    hData->SetMarkerStyle(20);
    hData->SetMarkerSize(1.15);
    hData->SetMarkerColor(kBlack);
    hData->SetLineColor(kBlack);

    hGauss->SetLineColor(kMagenta+1);
    hGauss->SetLineWidth(5);
    hGauss->SetLineStyle(1);
    hGauss->SetMarkerSize(0);

    hReso->SetLineColor(kBlue+1);
    hReso->SetLineWidth(6);
    hReso->SetLineStyle(2);
    hReso->SetMarkerSize(0);

    hRatioGauss->SetMarkerStyle(20);
    hRatioGauss->SetMarkerSize(1.0);
    hRatioGauss->SetMarkerColor(kMagenta+1);
    hRatioGauss->SetLineColor(kMagenta+1);

    hRatioReso->SetMarkerStyle(24);
    hRatioReso->SetMarkerSize(1.0);
    hRatioReso->SetMarkerColor(kBlue+1);
    hRatioReso->SetLineColor(kBlue+1);

    auto* c = new TCanvas(Form("cCATS3D_compare_%s", dir.Data()),
                          Form("CATS3D Gauss vs resonance %s", dir.Data()),
                          950, 900);

    auto* padTop = new TPad("padTop", "padTop", 0.0, 0.30, 1.0, 1.0);
    auto* padBot = new TPad("padBot", "padBot", 0.0, 0.0, 1.0, 0.30);
    padTop->SetBottomMargin(0.02);
    padTop->SetLeftMargin(0.15);
    padTop->SetRightMargin(0.04);
    padTop->SetTopMargin(0.04);
    padBot->SetTopMargin(0.03);
    padBot->SetBottomMargin(0.32);
    padBot->SetLeftMargin(0.15);
    padBot->SetRightMargin(0.04);
    padTop->Draw();
    padBot->Draw();

    padTop->cd();
    hData->SetTitle("");
    hData->GetXaxis()->SetLabelSize(0);
    hData->GetXaxis()->SetTitleSize(0);
    hData->GetYaxis()->SetTitle("C(q)");
    hData->GetYaxis()->SetTitleSize(0.060);
    hData->GetYaxis()->SetTitleOffset(1.00);
    hData->GetYaxis()->SetLabelSize(0.050);

    double yMin = hData->GetMinimum();
    double yMax = hData->GetMaximum();
    for(int i=1; i<=hGauss->GetNbinsX(); ++i){
        yMin = std::min(yMin, std::min(hGauss->GetBinContent(i),
                                       hReso->GetBinContent(i)));
        yMax = std::max(yMax, std::max(hGauss->GetBinContent(i),
                                       hReso->GetBinContent(i)));
    }
    hData->GetYaxis()->SetRangeUser(yMin - 0.08*(yMax-yMin),
                                    yMax + 0.24*(yMax-yMin));

    hData->Draw("E1");
    hReso->Draw("HIST SAME");
    hGauss->Draw("HIST SAME");
    hData->Draw("E1 SAME");

    auto* leg = new TLegend(0.34, 0.56, 0.93, 0.90);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.038);
    leg->AddEntry(hData, Form("data, %s slice", DirectionTitle(dir).Data()),
                  "lep");
    leg->AddEntry(hGauss,
                  "pure Gaussian, " + RadiusLabel(gaussSummary),
                  "l");
    leg->AddEntry(hReso,
                  Form("Gauss + resonance, #alpha=%.1f, %s",
                       resoSummary.Alpha, RadiusLabel(resoSummary).Data()),
                  "l");
    leg->Draw();

    padBot->cd();
    hRatioReso->SetTitle("");
    hRatioReso->GetXaxis()->SetTitle(AxisTitle(dir));
    hRatioReso->GetXaxis()->SetTitleSize(0.115);
    hRatioReso->GetXaxis()->SetTitleOffset(1.02);
    hRatioReso->GetXaxis()->SetLabelSize(0.095);
    hRatioReso->GetYaxis()->SetTitle("data / fit");
    hRatioReso->GetYaxis()->SetTitleSize(0.100);
    hRatioReso->GetYaxis()->SetTitleOffset(0.60);
    hRatioReso->GetYaxis()->SetLabelSize(0.085);
    hRatioReso->GetYaxis()->SetNdivisions(505);
    hRatioReso->GetYaxis()->SetRangeUser(0.90, 1.10);
    hRatioReso->Draw("E1");
    hRatioGauss->Draw("E1 SAME");

    auto* line = new TLine(hRatioReso->GetXaxis()->GetXmin(), 1.0,
                           hRatioReso->GetXaxis()->GetXmax(), 1.0);
    line->SetLineColor(kGray+2);
    line->SetLineStyle(2);
    line->Draw("SAME");

    c->cd();
    c->Update();

    TString saveName(outName);
    if(saveName.Length() == 0){
        saveName = Form("compare_CATS3D_gauss_vs_reso_%s.png",
                        dir.Data());
    }
    c->SaveAs(saveName);
}

void draw_compare_CATS3D_gauss_vs_reso_all()
{
    draw_compare_CATS3D_gauss_vs_reso(
        "out", true,
        "fit_CATS3D_gaussian_cent1_mt4_wtCoulomb.root",
        "fit_CATS3D_resonance_cent1_mt4_wtCoulomb_levyAlpha2.0.root",
        "compare_CATS3D_gauss_vs_reso_out.png");
    draw_compare_CATS3D_gauss_vs_reso(
        "side", true,
        "fit_CATS3D_gaussian_cent1_mt4_wtCoulomb.root",
        "fit_CATS3D_resonance_cent1_mt4_wtCoulomb_levyAlpha2.0.root",
        "compare_CATS3D_gauss_vs_reso_side.png");
    draw_compare_CATS3D_gauss_vs_reso(
        "long", true,
        "fit_CATS3D_gaussian_cent1_mt4_wtCoulomb.root",
        "fit_CATS3D_resonance_cent1_mt4_wtCoulomb_levyAlpha2.0.root",
        "compare_CATS3D_gauss_vs_reso_long.png");
}
