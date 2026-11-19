#include "TCanvas.h"
#include "TChain.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TNtuple.h"
#include "TStyle.h"

#include <cmath>
#include <iostream>

namespace {

double ReadFullValue(const int mt, const char* branch)
{
    TFile f(Form("fit_CATS3D_gaussian_cent1_mt%d_wtCoulomb.root", mt), "READ");
    auto* t = f.IsZombie() ? nullptr : dynamic_cast<TNtuple*>(f.Get("fit_summary"));
    if(!t || !t->GetLeaf(branch) || t->GetEntries()<1) return NAN;
    t->GetEntry(0);
    return t->GetLeaf(branch)->GetValue();
}

double BootstrapRms(const int mt, const char* branch)
{
    TChain t("fit_summary");
    for(int i=0; i<100; ++i){
        t.Add(Form("fit_CATS3D_gaussian_cent1_mt%d_boot_fullBaseline%d.root",
                   mt, i));
    }
    if(t.GetEntries()<2) return NAN;
    t.Draw(branch, "", "goff");
    return TMath::RMS(t.GetSelectedRows(), t.GetV1());
}

void StyleGraph(TGraphErrors& g, const Color_t color, const Style_t marker)
{
    g.SetLineColor(color);
    g.SetMarkerColor(color);
    g.SetMarkerStyle(marker);
    g.SetMarkerSize(1.35);
    g.SetLineWidth(2);
}

} // namespace

void draw_radii_vs_mt(const char* outName = "radii_vs_mt_cent1_stat.png")
{
    gStyle->SetOptStat(0);
    gStyle->SetEndErrorSize(5);

    constexpr int n = 4;
    const int mtIndex[n] = {2, 3, 4, 5};
    const double mt[n] = {0.45, 0.55, 0.65, 0.80};
    const double mtErr[n] = {0.05, 0.05, 0.05, 0.10};
    double rout[n], rside[n], rlong[n];
    double eout[n], eside[n], elong[n];

    for(int i=0; i<n; ++i){
        rout[i] = ReadFullValue(mtIndex[i], "Rout");
        rside[i] = ReadFullValue(mtIndex[i], "Rside");
        rlong[i] = ReadFullValue(mtIndex[i], "Rlong");
        eout[i] = BootstrapRms(mtIndex[i], "Rout");
        eside[i] = BootstrapRms(mtIndex[i], "Rside");
        elong[i] = BootstrapRms(mtIndex[i], "Rlong");
        std::cout << "mt" << mtIndex[i] << " " << mt[i]
                  << " Rout=" << rout[i] << "+-" << eout[i]
                  << " Rside=" << rside[i] << "+-" << eside[i]
                  << " Rlong=" << rlong[i] << "+-" << elong[i] << "\n";
    }

    TGraphErrors gOut(n, mt, rout, mtErr, eout);
    TGraphErrors gSide(n, mt, rside, mtErr, eside);
    TGraphErrors gLong(n, mt, rlong, mtErr, elong);
    StyleGraph(gOut, kRed+1, 20);
    StyleGraph(gSide, kBlue+1, 21);
    StyleGraph(gLong, kGreen+2, 33);

    TCanvas c("cRadiiVsMt", "Radii vs mT", 850, 700);
    c.SetLeftMargin(0.13);
    c.SetBottomMargin(0.13);
    c.SetRightMargin(0.04);
    c.SetTopMargin(0.05);
    gOut.SetTitle("");
    gOut.GetXaxis()->SetTitle("#it{m}_{T} (GeV/#it{c}^{2})");
    gOut.GetYaxis()->SetTitle("R (fm)");
    gOut.GetXaxis()->SetLimits(0.36, 0.93);
    gOut.GetYaxis()->SetRangeUser(3.1, 6.0);
    gOut.GetXaxis()->SetTitleSize(0.052);
    gOut.GetYaxis()->SetTitleSize(0.052);
    gOut.GetXaxis()->SetLabelSize(0.045);
    gOut.GetYaxis()->SetLabelSize(0.045);
    gOut.GetYaxis()->SetTitleOffset(1.05);
    gOut.Draw("APZ");
    gSide.Draw("PZ SAME");
    gLong.Draw("PZ SAME");

    TLegend leg(0.68, 0.72, 0.92, 0.92);
    leg.SetBorderSize(0);
    leg.SetFillStyle(0);
    leg.SetTextSize(0.043);
    leg.AddEntry(&gOut, "R_{out}", "lp");
    leg.AddEntry(&gSide, "R_{side}", "lp");
    leg.AddEntry(&gLong, "R_{long}", "lp");
    leg.Draw();

    c.SaveAs(outName);
    TString pdfName(outName);
    if(pdfName.EndsWith(".png")) pdfName.ReplaceAll(".png", ".pdf");
    c.SaveAs(pdfName);
}
