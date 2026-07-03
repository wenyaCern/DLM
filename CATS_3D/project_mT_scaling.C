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
#include <vector>

namespace {

double gNonFemtoQminMeV = 100.0;
double gQAxisScaleToMeV = 1000.0;

double ConstantNonFemtoTF3(double* x, double* par)
{
    if(std::fabs(x[0]) * gQAxisScaleToMeV < gNonFemtoQminMeV ||
       std::fabs(x[1]) * gQAxisScaleToMeV < gNonFemtoQminMeV ||
       std::fabs(x[2]) * gQAxisScaleToMeV < gNonFemtoQminMeV){
        TF1::RejectPoint();
    }
    return par[0];
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
    double qAxisScaleToMeV = 1000.0;
    double nonFemtoQminMeV = 100.0;
};

TF3* FitConstantNonFemto3D(TH3D* h, const MtProjConfig& cfg, const int mt)
{
    if(!h) return nullptr;

    int nBins = 0;

    for(int ix=1; ix<=h->GetNbinsX(); ++ix){
        for(int iy=1; iy<=h->GetNbinsY(); ++iy){
            for(int iz=1; iz<=h->GetNbinsZ(); ++iz){
                const double qOut = std::fabs(h->GetXaxis()->GetBinCenter(ix) *
                                             cfg.qAxisScaleToMeV);
                const double qSide = std::fabs(h->GetYaxis()->GetBinCenter(iy) *
                                              cfg.qAxisScaleToMeV);
                const double qLong = std::fabs(h->GetZaxis()->GetBinCenter(iz) *
                                              cfg.qAxisScaleToMeV);

                if(qOut  < cfg.nonFemtoQminMeV) continue;
                if(qSide < cfg.nonFemtoQminMeV) continue;
                if(qLong < cfg.nonFemtoQminMeV) continue;

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
    gQAxisScaleToMeV = cfg.qAxisScaleToMeV;

    TF3* fNorm = new TF3(Form("fConstNorm_mt%d", mt), ConstantNonFemtoTF3,
                         h->GetXaxis()->GetXmin(), h->GetXaxis()->GetXmax(),
                         h->GetYaxis()->GetXmin(), h->GetYaxis()->GetXmax(),
                         h->GetZaxis()->GetXmin(), h->GetZaxis()->GetXmax(),
                         1);
    fNorm->SetParName(0, "Norm");
    fNorm->SetParameter(0, 1.0);

    const int fitStatus = h->Fit(fNorm, "RQ0");

    std::cout << "mT bin " << mt
              << " constant non-femto TF3 fit norm = "
              << fNorm->GetParameter(0)
              << " +/- " << fNorm->GetParError(0)
              << " status=" << fitStatus
              << " from " << nBins
              << " bins with |q_out|, |q_side|, |q_long| > "
              << cfg.nonFemtoQminMeV << " MeV\n";

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

TH3D* BuildCk3D(TFile& fin, const MtProjConfig& cfg, const int mt, TF3*& fNorm)
{
    fNorm = nullptr;
    const TString sameName = Form("H3_same_cent%d_mt%d_qn%d_phi%d",
                                  cfg.cent, mt, cfg.qn, cfg.phi);
    const TString mixedName = Form("H3_mixed_cent%d_mt%d_phi%d",
                                   cfg.cent, mt, cfg.phi);

    TH3D* hSameIn = dynamic_cast<TH3D*>(fin.Get(sameName));
    TH3D* hMixedIn = dynamic_cast<TH3D*>(fin.Get(mixedName));
    if(!hSameIn || !hMixedIn){
        std::cout << "WARNING: missing " << sameName << " or " << mixedName << "\n";
        return nullptr;
    }

    std::unique_ptr<TH3D> hSame(dynamic_cast<TH3D*>(hSameIn->Clone(Form("same_mt%d_work", mt))));
    std::unique_ptr<TH3D> hMixed(dynamic_cast<TH3D*>(hMixedIn->Clone(Form("mixed_mt%d_work", mt))));
    hSame->SetDirectory(nullptr);
    hMixed->SetDirectory(nullptr);
    if(hSame->GetSumw2N() == 0) hSame->Sumw2();
    if(hMixed->GetSumw2N() == 0) hMixed->Sumw2();

    if(cfg.rebin > 1){
        hSame->Rebin3D(cfg.rebin, cfg.rebin, cfg.rebin);
        hMixed->Rebin3D(cfg.rebin, cfg.rebin, cfg.rebin);
    }

    TH3D* hCk = dynamic_cast<TH3D*>(hSame->Clone(Form("hCk3D_mt%d", mt)));
    hCk->SetDirectory(nullptr);
    hCk->Divide(hSame.get(), hMixed.get(), 1.0, 1.0);
    fNorm = FitConstantNonFemto3D(hCk, cfg, mt);
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

    std::vector<TH1D*> projOut;
    std::vector<TH1D*> projSide;
    std::vector<TH1D*> projLong;
    std::vector<int> goodMtBins;

    for(const int mt : cfg.mtBins){
        TF3* fNormRaw = nullptr;
        std::unique_ptr<TH3D> hCk(BuildCk3D(fin, cfg, mt, fNormRaw));
        std::unique_ptr<TF3> fNorm(fNormRaw);
        if(!hCk) continue;

        fout.cd();
        if(fNorm) fNorm->Write("", TObject::kOverwrite);
        hCk->Write("", TObject::kOverwrite);

        TH1D* hOut = ProjectFemtoSlice(hCk.get(), "out",
                                       Form("hCk_proj_out_mt%d", mt),
                                       cfg.sliceCut);
        TH1D* hSide = ProjectFemtoSlice(hCk.get(), "side",
                                        Form("hCk_proj_side_mt%d", mt),
                                        cfg.sliceCut);
        TH1D* hLong = ProjectFemtoSlice(hCk.get(), "long",
                                        Form("hCk_proj_long_mt%d", mt),
                                        cfg.sliceCut);

        if(hOut) hOut->Write("", TObject::kOverwrite);
        if(hSide) hSide->Write("", TObject::kOverwrite);
        if(hLong) hLong->Write("", TObject::kOverwrite);

        projOut.push_back(hOut);
        projSide.push_back(hSide);
        projLong.push_back(hLong);
        goodMtBins.push_back(mt);
    }

    fout.Close();

    // DrawOverlay(projOut, goodMtBins, "out", cfg.outputPrefix + "_out.png");
    // DrawOverlay(projSide, goodMtBins, "side", cfg.outputPrefix + "_side.png");
    // DrawOverlay(projLong, goodMtBins, "long", cfg.outputPrefix + "_long.png");

    std::cout << "Wrote " << cfg.outputFile << "\n";
    std::cout << "Wrote " << cfg.outputPrefix << "_out/side/long.png\n";
}
