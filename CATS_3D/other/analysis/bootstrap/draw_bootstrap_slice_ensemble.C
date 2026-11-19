#include "TCanvas.h"
#include "TFile.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TString.h"

#include <memory>
#include <vector>

namespace {

const char* AxisTitle(const char* direction)
{
    const TString d(direction);
    if(d=="out") return "q_{out} (GeV/c)";
    if(d=="side") return "q_{side} (GeV/c)";
    return "q_{long} (GeV/c)";
}

std::vector<std::unique_ptr<TH1D>> LoadEnsemble(const char* pattern,
                                                const char* histPattern,
                                                const char* direction,
                                                const int nSamples)
{
    std::vector<std::unique_ptr<TH1D>> result;
    for(int i=0; i<nSamples; ++i){
        TFile f(Form(pattern, i), "READ");
        auto* input = f.IsZombie() ? nullptr
                      : dynamic_cast<TH1D*>(f.Get(Form(histPattern, direction)));
        if(!input) continue;
        auto* copy = dynamic_cast<TH1D*>(input->Clone(
            Form("ensemble_%s_%d", direction, i)));
        copy->SetDirectory(nullptr);
        result.emplace_back(copy);
    }
    return result;
}

void DrawDirection(TPad* pad,
                   const char* fullFileName,
                   const char* bootstrapPattern,
                   const char* direction,
                   const char* ensembleHistPattern,
                   const char* fullHistPattern,
                   const int nSamples,
                   const bool drawData)
{
    pad->cd();
    pad->SetLeftMargin(0.14);
    pad->SetRightMargin(0.025);
    pad->SetBottomMargin(0.14);
    pad->SetTopMargin(0.06);

    auto ensemble = LoadEnsemble(bootstrapPattern, ensembleHistPattern,
                                 direction, nSamples);
    if(ensemble.empty()) return;

    TFile fullFile(fullFileName, "READ");
    auto* fullInput = dynamic_cast<TH1D*>(
        fullFile.Get(Form(fullHistPattern, direction)));
    if(!fullInput) return;
    std::unique_ptr<TH1D> full(dynamic_cast<TH1D*>(fullInput->Clone(
        Form("full_%s", direction))));
    full->SetDirectory(nullptr);

    TH1D* frame = ensemble.front().get();
    frame->SetTitle("");
    frame->GetXaxis()->SetTitle(AxisTitle(direction));
    frame->GetYaxis()->SetTitle("C(q)");
    frame->GetXaxis()->SetRangeUser(-0.2, 0.2);
    frame->GetYaxis()->SetRangeUser(0.98, 1.35);
    frame->GetXaxis()->SetTitleSize(0.052);
    frame->GetYaxis()->SetTitleSize(0.052);
    frame->GetXaxis()->SetLabelSize(0.044);
    frame->GetYaxis()->SetLabelSize(0.044);
    frame->GetYaxis()->SetTitleOffset(1.15);
    frame->SetLineColorAlpha(drawData ? kGray+2 : kAzure+2, 0.13);
    frame->SetLineWidth(1);
    frame->SetMarkerSize(0);
    auto* frameDraw = dynamic_cast<TH1D*>(frame->DrawCopy("HIST"));
    for(size_t i=1; i<ensemble.size(); ++i){
        ensemble[i]->SetLineColorAlpha(drawData ? kGray+2 : kAzure+2, 0.13);
        ensemble[i]->SetLineWidth(1);
        ensemble[i]->SetMarkerSize(0);
        ensemble[i]->DrawCopy("HIST SAME");
    }

    TH1D* fullDraw = nullptr;
    if(drawData){
        full->SetMarkerStyle(20);
        full->SetMarkerSize(0.85);
        full->SetMarkerColor(kBlack);
        full->SetLineColor(kBlack);
        fullDraw = dynamic_cast<TH1D*>(full->DrawCopy("E1 SAME"));
    }
    else {
        full->SetLineColor(kBlue+2);
        full->SetLineWidth(3);
        full->SetMarkerSize(0);
        fullDraw = dynamic_cast<TH1D*>(full->DrawCopy("HIST SAME"));
    }

    TLatex label;
    label.SetNDC();
    label.SetTextSize(0.050);
    label.DrawLatex(0.18, 0.88, Form("q_{%s}", direction));

    auto* legend = new TLegend(0.50, 0.76, 0.94, 0.91);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->SetTextSize(0.037);
    legend->AddEntry(frameDraw, Form("%zu bootstrap replicas", ensemble.size()), "l");
    legend->AddEntry(fullDraw, drawData ? "full data" : "full-data best fit",
                     drawData ? "lep" : "l");
    legend->Draw();
}

void DrawEnsembleCanvas(const char* fullFileName,
                        const char* bootstrapPattern,
                        const char* outName,
                        const int nSamples,
                        const bool drawData)
{
    TCanvas canvas(drawData ? "cBootstrapData" : "cBootstrapFits",
                   "bootstrap ensemble", 1500, 520);
    canvas.Divide(3, 1, 0.001, 0.001);
    const char* directions[3] = {"out", "side", "long"};
    for(int i=0; i<3; ++i){
        DrawDirection(dynamic_cast<TPad*>(canvas.cd(i+1)), fullFileName,
                      bootstrapPattern, directions[i],
                      drawData ? "hDataCk3D_proj_%s" : "hCATS3D_best_proj_%s",
                      drawData ? "hDataCk3D_proj_%s" : "hCATS3D_best_proj_%s",
                      nSamples, drawData);
    }
    canvas.SaveAs(outName);
    TString pdf(outName);
    if(pdf.EndsWith(".png")) pdf.ReplaceAll(".png", ".pdf");
    canvas.SaveAs(pdf);
}

} // namespace

void draw_bootstrap_slice_ensemble(
    const char* fullFileName = "fit_CATS3D_gaussian_cent1_mt2_wtCoulomb.root",
    const char* bootstrapPattern =
        "fit_CATS3D_gaussian_cent1_mt2_boot_fullBaseline%d.root",
    const int nSamples = 100)
{
    gStyle->SetOptStat(0);
    DrawEnsembleCanvas(fullFileName, bootstrapPattern,
                       "mt2_bootstrap_100_data_slices.png", nSamples, true);
    DrawEnsembleCanvas(fullFileName, bootstrapPattern,
                       "mt2_bootstrap_100_fit_lines.png", nSamples, false);
}
