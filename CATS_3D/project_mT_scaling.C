#include "TCanvas.h"
#include "TFile.h"
#include "TF1.h"
#include "TF3.h"
#include "TH1D.h"
#include "TH3D.h"
#include "TLegend.h"
#include "TMath.h"
#include "TString.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <memory>
#include <utility>
#include <vector>

namespace {

double gNonFemtoQminMeV = 100.0;
double gNonFemtoQmaxMeV = 200.0;
double gQAxisScaleToMeV = 1000.0;
int gNormMode = 0;

double ConstantNonFemtoTF3(double* x, double* par)
{
    const double qOut = std::fabs(x[0]) * gQAxisScaleToMeV;
    const double qSide = std::fabs(x[1]) * gQAxisScaleToMeV;
    const double qLong = std::fabs(x[2]) * gQAxisScaleToMeV;
    const double kStar = 0.5*std::sqrt(qOut*qOut + qSide*qSide + qLong*qLong);

    bool accept = false;
    if(gNormMode == 0){
        accept = qOut >= gNonFemtoQminMeV &&
                 qSide >= gNonFemtoQminMeV &&
                 qLong >= gNonFemtoQminMeV;
    }
    else if(gNormMode == 1){
        accept = kStar >= gNonFemtoQminMeV &&
                 kStar <= gNonFemtoQmaxMeV;
    }
    else if(gNormMode == 2){
        accept = qOut >= gNonFemtoQminMeV ||
                 qSide >= gNonFemtoQminMeV ||
                 qLong >= gNonFemtoQminMeV;
    }
    else {
        accept = true;
    }

    if(!accept){
        TF1::RejectPoint();
    }
    return par[0];
}

const char* NormModeName(const int mode)
{
    if(mode == 0) return "corner";
    if(mode == 1) return "kstar";
    if(mode == 2) return "outerBox";
    if(mode == 3) return "none";
    return "unknown";
}

struct MtProjConfig {
    TString inputFile = "File/out_3DCF_mT_cent_qn_phi_tofOnly.root";
    TString outputFile = "OutputFiles/mT_scaling_projection_norm.root";
    TString outputPrefix = "OutputFiles/mT_scaling";

    int cent = 0;
    int qn = 0;
    int phi = 0;
    std::vector<int> mtBins = {0, 1, 2, 3, 4, 5, 6, 7};

    int rebin = 1;
    double sliceCut = 0.04;
    std::vector<double> sliceCutScan = {0.02, 0.04, 0.06, 0.08};
    double qAxisScaleToMeV = 1000.0;
    double nonFemtoQminMeV = 100.0;
    double nonFemtoQmaxMeV = 200.0;
    std::vector<int> normModes = {0, 1, 2, 3};
};

bool IsNormBin(const TH3D* h,
               const int ix,
               const int iy,
               const int iz,
               const MtProjConfig& cfg,
               const int normMode)
{
    const double qOut = std::fabs(h->GetXaxis()->GetBinCenter(ix) *
                                  cfg.qAxisScaleToMeV);
    const double qSide = std::fabs(h->GetYaxis()->GetBinCenter(iy) *
                                   cfg.qAxisScaleToMeV);
    const double qLong = std::fabs(h->GetZaxis()->GetBinCenter(iz) *
                                   cfg.qAxisScaleToMeV);
    const double kStar = 0.5*std::sqrt(qOut*qOut + qSide*qSide + qLong*qLong);

    if(normMode == 0){
        return qOut >= cfg.nonFemtoQminMeV &&
               qSide >= cfg.nonFemtoQminMeV &&
               qLong >= cfg.nonFemtoQminMeV;
    }
    if(normMode == 1){
        return kStar >= cfg.nonFemtoQminMeV &&
               kStar <= cfg.nonFemtoQmaxMeV;
    }
    if(normMode == 2){
        return qOut >= cfg.nonFemtoQminMeV ||
               qSide >= cfg.nonFemtoQminMeV ||
               qLong >= cfg.nonFemtoQminMeV;
    }
    return false;
}

TF3* FitConstantNonFemto3D(TH3D* h, const MtProjConfig& cfg, const int mt,
                           const int normMode)
{
    if(!h) return nullptr;
    if(normMode == 3) return nullptr;

    int nBins = 0;

    for(int ix=1; ix<=h->GetNbinsX(); ++ix){
        for(int iy=1; iy<=h->GetNbinsY(); ++iy){
            for(int iz=1; iz<=h->GetNbinsZ(); ++iz){
                if(!IsNormBin(h, ix, iy, iz, cfg, normMode)) continue;

                const double c = h->GetBinContent(ix, iy, iz);
                const double e = h->GetBinError(ix, iy, iz);
                if(!std::isfinite(c) || c <= 0.0) continue;
                ++nBins;
            }
        }
    }

    if(nBins <= 0){
        std::cout << "WARNING: no non-femto normalization bins for "
                  << h->GetName() << "\n";
        return nullptr;
    }

    gNonFemtoQminMeV = cfg.nonFemtoQminMeV;
    gNonFemtoQmaxMeV = cfg.nonFemtoQmaxMeV;
    gQAxisScaleToMeV = cfg.qAxisScaleToMeV;
    gNormMode = normMode;

    TF3* fNorm = new TF3(Form("fConstNorm_%s_mt%d", NormModeName(normMode), mt),
                         ConstantNonFemtoTF3,
                         h->GetXaxis()->GetXmin(), h->GetXaxis()->GetXmax(),
                         h->GetYaxis()->GetXmin(), h->GetYaxis()->GetXmax(),
                         h->GetZaxis()->GetXmin(), h->GetZaxis()->GetXmax(),
                         1);
    fNorm->SetParName(0, "Norm");
    fNorm->SetParameter(0, 1.0);

    const int fitStatus = h->Fit(fNorm, "RQ0");

    std::cout << "mT bin " << mt
              << " normMode=" << NormModeName(normMode)
              << " constant non-femto TF3 fit norm = "
              << fNorm->GetParameter(0)
              << " +/- " << fNorm->GetParError(0)
              << " status=" << fitStatus
              << " from " << nBins
              << " bins\n";

    return fNorm;
}

void DivideByConstantNorm(TH3D* h, const TF3* fNorm)
{
    if(!h || !fNorm) return;
    const double norm = fNorm->GetParameter(0);
    if(norm > 0.0 && std::isfinite(norm)){
        h->Scale(1.0 / norm);
    }
}

bool LoadSameMixed3D(TFile& fin,
                     const MtProjConfig& cfg,
                     const int mt,
                     const TString& tag,
                     std::unique_ptr<TH3D>& hSame,
                     std::unique_ptr<TH3D>& hMixed)
{
    const TString sameName = Form("H3_same_cent%d_mt%d_qn%d_phi%d",
                                  cfg.cent, mt, cfg.qn, cfg.phi);
    const TString mixedName = Form("H3_mixed_cent%d_mt%d_phi%d",
                                   cfg.cent, mt, cfg.phi);

    TH3D* hSameIn = dynamic_cast<TH3D*>(fin.Get(sameName));
    TH3D* hMixedIn = dynamic_cast<TH3D*>(fin.Get(mixedName));
    if(!hSameIn || !hMixedIn){
        std::cout << "WARNING: missing " << sameName << " or " << mixedName << "\n";
        return false;
    }

    hSame.reset(dynamic_cast<TH3D*>(hSameIn->Clone(Form("same_%s_mt%d_work", tag.Data(), mt))));
    hMixed.reset(dynamic_cast<TH3D*>(hMixedIn->Clone(Form("mixed_%s_mt%d_work", tag.Data(), mt))));
    hSame->SetDirectory(nullptr);
    hMixed->SetDirectory(nullptr);
    if(hSame->GetSumw2N() == 0) hSame->Sumw2();
    if(hMixed->GetSumw2N() == 0) hMixed->Sumw2();

    if(cfg.rebin > 1){
        hSame->Rebin3D(cfg.rebin, cfg.rebin, cfg.rebin);
        hMixed->Rebin3D(cfg.rebin, cfg.rebin, cfg.rebin);
    }

    return true;
}

TH3D* BuildCk3D(TH3D* hSame, TH3D* hMixed, const MtProjConfig& cfg,
                const int mt, const int normMode, TF3*& fNorm)
{
    fNorm = nullptr;
    if(!hSame || !hMixed) return nullptr;

    TH3D* hCk = dynamic_cast<TH3D*>(hSame->Clone(Form("hCk3D_%s_mt%d", NormModeName(normMode), mt)));
    hCk->SetDirectory(nullptr);
    hCk->Divide(hSame, hMixed, 1.0, 1.0);
    fNorm = FitConstantNonFemto3D(hCk, cfg, mt, normMode);
    DivideByConstantNorm(hCk, fNorm);
    return hCk;
}

TH1D* ProjectFemtoSlice(TH3D* h3, const TString& axis, const TString& name,
                        const double cutAbs)
{
    if(!h3) return nullptr;
    TH1D* h1 = nullptr;

    if(axis == "out"){
        const int y1 = h3->GetYaxis()->FindBin(-cutAbs);
        const int y2 = h3->GetYaxis()->FindBin( cutAbs);
        const int z1 = h3->GetZaxis()->FindBin(-cutAbs);
        const int z2 = h3->GetZaxis()->FindBin( cutAbs);
        h1 = dynamic_cast<TH1D*>(h3->ProjectionX(name, y1, y2, z1, z2, "e"));
        h1->Scale(1.0 / double((y2-y1+1)*(z2-z1+1)));
    }
    else if(axis == "side"){
        const int x1 = h3->GetXaxis()->FindBin(-cutAbs);
        const int x2 = h3->GetXaxis()->FindBin( cutAbs);
        const int z1 = h3->GetZaxis()->FindBin(-cutAbs);
        const int z2 = h3->GetZaxis()->FindBin( cutAbs);
        h1 = dynamic_cast<TH1D*>(h3->ProjectionY(name, x1, x2, z1, z2, "e"));
        h1->Scale(1.0 / double((x2-x1+1)*(z2-z1+1)));
    }
    else if(axis == "long"){
        const int x1 = h3->GetXaxis()->FindBin(-cutAbs);
        const int x2 = h3->GetXaxis()->FindBin( cutAbs);
        const int y1 = h3->GetYaxis()->FindBin(-cutAbs);
        const int y2 = h3->GetYaxis()->FindBin( cutAbs);
        h1 = dynamic_cast<TH1D*>(h3->ProjectionZ(name, x1, x2, y1, y2, "e"));
        h1->Scale(1.0 / double((x2-x1+1)*(y2-y1+1)));
    }

    if(h1){
        h1->SetDirectory(nullptr);
        h1->SetTitle(Form("%s;#it{q}_{%s} (GeV/#it{c});C(#it{q})",
                          name.Data(), axis.Data()));
    }

    TH1D* h1_copy = new TH1D(Form("scaled_%s;#it{q}_{%s} (GeV/#it{c});C(#it{q})",
                          name.Data(), axis.Data()),Form("scaled_%s;#it{q}_{%s} (GeV/#it{c});C(#it{q})",
                          name.Data(), axis.Data()), h1->GetNbinsX(), -200, 200);
    for (int i(0); i<h1->GetNbinsX(); ++i){
        h1_copy->SetBinContent(i+1, h1->GetBinContent(i+1));
        h1_copy->SetBinError(i+1, h1->GetBinError(i+1));
    }
    return h1_copy;
}

TH1D* ProjectSameMixedRatioSlice(TH3D* hSame,
                                 TH3D* hMixed,
                                 const TString& axis,
                                 const TString& name,
                                 const double cutAbs,
                                 const TF3* fNorm)
{
    if(!hSame || !hMixed) return nullptr;

    TH1D* hNum = nullptr;
    TH1D* hDen = nullptr;

    if(axis == "out"){
        const int y1 = hSame->GetYaxis()->FindBin(-cutAbs);
        const int y2 = hSame->GetYaxis()->FindBin( cutAbs);
        const int z1 = hSame->GetZaxis()->FindBin(-cutAbs);
        const int z2 = hSame->GetZaxis()->FindBin( cutAbs);
        hNum = dynamic_cast<TH1D*>(hSame->ProjectionX(name + "_num", y1, y2, z1, z2, "e"));
        hDen = dynamic_cast<TH1D*>(hMixed->ProjectionX(name + "_den", y1, y2, z1, z2, "e"));
    }
    else if(axis == "side"){
        const int x1 = hSame->GetXaxis()->FindBin(-cutAbs);
        const int x2 = hSame->GetXaxis()->FindBin( cutAbs);
        const int z1 = hSame->GetZaxis()->FindBin(-cutAbs);
        const int z2 = hSame->GetZaxis()->FindBin( cutAbs);
        hNum = dynamic_cast<TH1D*>(hSame->ProjectionY(name + "_num", x1, x2, z1, z2, "e"));
        hDen = dynamic_cast<TH1D*>(hMixed->ProjectionY(name + "_den", x1, x2, z1, z2, "e"));
    }
    else if(axis == "long"){
        const int x1 = hSame->GetXaxis()->FindBin(-cutAbs);
        const int x2 = hSame->GetXaxis()->FindBin( cutAbs);
        const int y1 = hSame->GetYaxis()->FindBin(-cutAbs);
        const int y2 = hSame->GetYaxis()->FindBin( cutAbs);
        hNum = dynamic_cast<TH1D*>(hSame->ProjectionZ(name + "_num", x1, x2, y1, y2, "e"));
        hDen = dynamic_cast<TH1D*>(hMixed->ProjectionZ(name + "_den", x1, x2, y1, y2, "e"));
    }

    if(!hNum || !hDen) return nullptr;
    hNum->SetDirectory(nullptr);
    hDen->SetDirectory(nullptr);

    TH1D* hRatio = dynamic_cast<TH1D*>(hNum->Clone(name));
    hRatio->SetDirectory(nullptr);
    hRatio->Divide(hNum, hDen, 1.0, 1.0);

    if(fNorm){
        const double norm = fNorm->GetParameter(0);
        if(norm > 0.0 && std::isfinite(norm)) hRatio->Scale(1.0 / norm);
    }

    hRatio->SetTitle(Form("%s;#it{q}_{%s} (GeV/#it{c});C(#it{q})",
                          name.Data(), axis.Data()));

    TH1D* hRatioCopy = new TH1D(Form("scaled_%s;#it{q}_{%s} (GeV/#it{c});C(#it{q})",
                                     name.Data(), axis.Data()),
                                Form("scaled_%s;#it{q}_{%s} (GeV/#it{c});C(#it{q})",
                                     name.Data(), axis.Data()),
                                hRatio->GetNbinsX(), -200, 200);
    for(int i=1; i<=hRatio->GetNbinsX(); ++i){
        hRatioCopy->SetBinContent(i, hRatio->GetBinContent(i));
        hRatioCopy->SetBinError(i, hRatio->GetBinError(i));
    }

    delete hNum;
    delete hDen;
    delete hRatio;
    return hRatioCopy;
}

double PeakNearZero(const TH1D* h)
{
    if(!h) return 0.0;
    double sumW = 0.0;
    double sumY = 0.0;
    for(int i=1; i<=h->GetNbinsX(); ++i){
        const double x = h->GetBinCenter(i);
        if(std::fabs(x) > 20.0) continue;
        const double y = h->GetBinContent(i);
        const double e = h->GetBinError(i);
        if(!std::isfinite(y) || y <= 0.0) continue;
        const double w = (e > 0.0 && std::isfinite(e)) ? 1.0/(e*e) : 1.0;
        sumW += w;
        sumY += w*y;
    }
    return sumW > 0.0 ? sumY/sumW : 0.0;
}

void DrawOverlay(const std::vector<TH1D*>& histos,
                 const std::vector<int>& mtBins,
                 const TString& axis,
                 const TString& pngName)
{
    if(histos.empty()) return;

    const int colors[] = {kBlack, kRed+1, kBlue+1, kGreen+2, kMagenta+1,
                          kOrange+7, kCyan+2, kViolet+1};

    TCanvas c(Form("c_%s", axis.Data()), Form("c_%s", axis.Data()), 900, 650);
    c.SetTicks(1, 1);

    TLegend leg(0.62, 0.70, 0.88, 0.88);
    leg.SetBorderSize(0);
    leg.SetFillStyle(0);

    double yMin = 1e9;
    double yMax = -1e9;
    for(TH1D* h : histos){
        if(!h) continue;
        yMin = std::min(yMin, h->GetMinimum());
        yMax = std::max(yMax, h->GetMaximum());
    }
    if(yMin >= yMax){
        yMin = 0.8;
        yMax = 1.2;
    }

    for(size_t i=0; i<histos.size(); ++i){
        TH1D* h = histos[i];
        if(!h) continue;
        h->SetMarkerStyle(20 + int(i % 7));
        h->SetMarkerSize(0.8);
        h->SetMarkerColor(colors[i % 8]);
        h->SetLineColor(colors[i % 8]);
        h->SetLineWidth(2);
        h->SetMinimum(yMin - 0.15*(yMax-yMin));
        h->SetMaximum(yMax + 0.25*(yMax-yMin));

        h->Draw(i == 0 ? "e1" : "e1 same");
        leg.AddEntry(h, Form("mT bin %d", mtBins[i]), "lep");
    }

    leg.Draw();
    c.SaveAs(pngName);
}

void DrawCutScanOverlay(const std::vector<TH1D*>& histos,
                        const std::vector<double>& cuts,
                        const int mt,
                        const TString& axis,
                        const TString& pngName)
{
    if(histos.empty()) return;

    const int colors[] = {kBlack, kRed+1, kBlue+1, kGreen+2, kMagenta+1,
                          kOrange+7, kCyan+2, kViolet+1};

    TCanvas c(Form("c_cutscan_%s_mt%d", axis.Data(), mt),
              Form("c_cutscan_%s_mt%d", axis.Data(), mt), 900, 650);
    c.SetTicks(1, 1);

    TLegend leg(0.58, 0.68, 0.88, 0.88);
    leg.SetBorderSize(0);
    leg.SetFillStyle(0);

    double yMin = 1e9;
    double yMax = -1e9;
    for(TH1D* h : histos){
        if(!h) continue;
        yMin = std::min(yMin, h->GetMinimum());
        yMax = std::max(yMax, h->GetMaximum());
    }
    if(yMin >= yMax){
        yMin = 0.8;
        yMax = 1.2;
    }

    for(size_t i=0; i<histos.size(); ++i){
        TH1D* h = histos[i];
        if(!h) continue;
        h->SetMarkerStyle(20 + int(i % 7));
        h->SetMarkerSize(0.8);
        h->SetMarkerColor(colors[i % 8]);
        h->SetLineColor(colors[i % 8]);
        h->SetLineWidth(2);
        h->SetMinimum(yMin - 0.15*(yMax-yMin));
        h->SetMaximum(yMax + 0.25*(yMax-yMin));
        h->Draw(i == 0 ? "e1" : "e1 same");
        leg.AddEntry(h, Form("|q_{other}| < %.0f MeV", cuts[i]*1000.0), "lep");
    }

    leg.Draw();
    c.SaveAs(pngName);
}

} // namespace

void project_mT_scaling()
{
    MtProjConfig cfg;

    TFile fin(cfg.inputFile, "READ");
    if(fin.IsZombie()){
        std::cout << "ERROR: cannot open " << cfg.inputFile << "\n";
        return;
    }

    TFile fout(cfg.outputFile, "RECREATE");

    for(const int normMode : cfg.normModes){
        std::vector<TH1D*> projOut;
        std::vector<TH1D*> projSide;
        std::vector<TH1D*> projLong;
        std::vector<int> goodMtBins;

        std::cout << "\n=== normalization mode: "
                  << NormModeName(normMode) << " ===\n";

        for(const int mt : cfg.mtBins){
            std::unique_ptr<TH3D> hSame;
            std::unique_ptr<TH3D> hMixed;
            if(!LoadSameMixed3D(fin, cfg, mt, NormModeName(normMode), hSame, hMixed)){
                continue;
            }

            TF3* fNormRaw = nullptr;
            std::unique_ptr<TH3D> hCk(BuildCk3D(hSame.get(), hMixed.get(), cfg,
                                                mt, normMode, fNormRaw));
            std::unique_ptr<TF3> fNorm(fNormRaw);
            if(!hCk) continue;

            fout.cd();
            if(fNorm) fNorm->Write("", TObject::kOverwrite);
            hCk->Write("", TObject::kOverwrite);

            TH1D* hOut = ProjectFemtoSlice(
                hCk.get(), "out",
                Form("hCk_proj_out_%s_mt%d", NormModeName(normMode), mt),
                cfg.sliceCut);
            TH1D* hSide = ProjectFemtoSlice(
                hCk.get(), "side",
                Form("hCk_proj_side_%s_mt%d", NormModeName(normMode), mt),
                cfg.sliceCut);
            TH1D* hLong = ProjectFemtoSlice(
                hCk.get(), "long",
                Form("hCk_proj_long_%s_mt%d", NormModeName(normMode), mt),
                cfg.sliceCut);

            if(hOut) hOut->Write("", TObject::kOverwrite);
            if(hSide) hSide->Write("", TObject::kOverwrite);
            if(hLong) hLong->Write("", TObject::kOverwrite);

            TH1D* hOutRatio = ProjectSameMixedRatioSlice(
                hSame.get(), hMixed.get(), "out",
                Form("hCk_ratioProj_out_%s_mt%d", NormModeName(normMode), mt),
                cfg.sliceCut, fNorm.get());
            TH1D* hSideRatio = ProjectSameMixedRatioSlice(
                hSame.get(), hMixed.get(), "side",
                Form("hCk_ratioProj_side_%s_mt%d", NormModeName(normMode), mt),
                cfg.sliceCut, fNorm.get());
            TH1D* hLongRatio = ProjectSameMixedRatioSlice(
                hSame.get(), hMixed.get(), "long",
                Form("hCk_ratioProj_long_%s_mt%d", NormModeName(normMode), mt),
                cfg.sliceCut, fNorm.get());

            if(hOutRatio) hOutRatio->Write("", TObject::kOverwrite);
            if(hSideRatio) hSideRatio->Write("", TObject::kOverwrite);
            if(hLongRatio) hLongRatio->Write("", TObject::kOverwrite);

            std::cout << "projection order check mode=" << NormModeName(normMode)
                      << " mt=" << mt
                      << " cut=" << cfg.sliceCut
                      << " avgPeak(out,side,long)=("
                      << PeakNearZero(hOut) << ", "
                      << PeakNearZero(hSide) << ", "
                      << PeakNearZero(hLong) << ")"
                      << " ratioPeak(out,side,long)=("
                      << PeakNearZero(hOutRatio) << ", "
                      << PeakNearZero(hSideRatio) << ", "
                      << PeakNearZero(hLongRatio) << ")\n";

            projOut.push_back(hOut);
            projSide.push_back(hSide);
            projLong.push_back(hLong);
            goodMtBins.push_back(mt);

            std::vector<TH1D*> cutScanOut;
            std::vector<TH1D*> cutScanSide;
            std::vector<TH1D*> cutScanLong;

            for(const double cut : cfg.sliceCutScan){
                TH1D* hOutCut = ProjectFemtoSlice(
                    hCk.get(), "out",
                    Form("hCk_proj_out_%s_mt%d_cut%.0fMeV",
                         NormModeName(normMode), mt, cut*1000.0),
                    cut);
                TH1D* hSideCut = ProjectFemtoSlice(
                    hCk.get(), "side",
                    Form("hCk_proj_side_%s_mt%d_cut%.0fMeV",
                         NormModeName(normMode), mt, cut*1000.0),
                    cut);
                TH1D* hLongCut = ProjectFemtoSlice(
                    hCk.get(), "long",
                    Form("hCk_proj_long_%s_mt%d_cut%.0fMeV",
                         NormModeName(normMode), mt, cut*1000.0),
                    cut);
                TH1D* hOutRatioCut = ProjectSameMixedRatioSlice(
                    hSame.get(), hMixed.get(), "out",
                    Form("hCk_ratioProj_out_%s_mt%d_cut%.0fMeV",
                         NormModeName(normMode), mt, cut*1000.0),
                    cut, fNorm.get());
                TH1D* hSideRatioCut = ProjectSameMixedRatioSlice(
                    hSame.get(), hMixed.get(), "side",
                    Form("hCk_ratioProj_side_%s_mt%d_cut%.0fMeV",
                         NormModeName(normMode), mt, cut*1000.0),
                    cut, fNorm.get());
                TH1D* hLongRatioCut = ProjectSameMixedRatioSlice(
                    hSame.get(), hMixed.get(), "long",
                    Form("hCk_ratioProj_long_%s_mt%d_cut%.0fMeV",
                         NormModeName(normMode), mt, cut*1000.0),
                    cut, fNorm.get());

                if(hOutCut) hOutCut->Write("", TObject::kOverwrite);
                if(hSideCut) hSideCut->Write("", TObject::kOverwrite);
                if(hLongCut) hLongCut->Write("", TObject::kOverwrite);
                if(hOutRatioCut) hOutRatioCut->Write("", TObject::kOverwrite);
                if(hSideRatioCut) hSideRatioCut->Write("", TObject::kOverwrite);
                if(hLongRatioCut) hLongRatioCut->Write("", TObject::kOverwrite);

                std::cout << "norm scan mode=" << NormModeName(normMode)
                          << " mt=" << mt
                          << " cut=" << cut
                          << " avgPeak out=" << PeakNearZero(hOutCut)
                          << " side=" << PeakNearZero(hSideCut)
                          << " long=" << PeakNearZero(hLongCut)
                          << " ratioPeak out=" << PeakNearZero(hOutRatioCut)
                          << " side=" << PeakNearZero(hSideRatioCut)
                          << " long=" << PeakNearZero(hLongRatioCut)
                          << "\n";

                cutScanOut.push_back(hOutCut);
                cutScanSide.push_back(hSideCut);
                cutScanLong.push_back(hLongCut);
            }

            DrawCutScanOverlay(
                cutScanOut, cfg.sliceCutScan, mt, "out",
                Form("%s_%s_cutScan_out_mt%d.png",
                     cfg.outputPrefix.Data(), NormModeName(normMode), mt));
            DrawCutScanOverlay(
                cutScanSide, cfg.sliceCutScan, mt, "side",
                Form("%s_%s_cutScan_side_mt%d.png",
                     cfg.outputPrefix.Data(), NormModeName(normMode), mt));
            DrawCutScanOverlay(
                cutScanLong, cfg.sliceCutScan, mt, "long",
                Form("%s_%s_cutScan_long_mt%d.png",
                     cfg.outputPrefix.Data(), NormModeName(normMode), mt));
        }

        DrawOverlay(projOut, goodMtBins, "out",
                    Form("%s_%s_out.png", cfg.outputPrefix.Data(),
                         NormModeName(normMode)));
        DrawOverlay(projSide, goodMtBins, "side",
                    Form("%s_%s_side.png", cfg.outputPrefix.Data(),
                         NormModeName(normMode)));
        DrawOverlay(projLong, goodMtBins, "long",
                    Form("%s_%s_long.png", cfg.outputPrefix.Data(),
                         NormModeName(normMode)));
    }

    fout.Close();

    std::cout << "Wrote " << cfg.outputFile << "\n";
    std::cout << "Wrote normalization comparison png files under "
              << cfg.outputPrefix << "_*.png\n";
}
