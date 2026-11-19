#include "TCanvas.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH3D.h"
#include "TLegend.h"
#include "TString.h"
#include "TStyle.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <memory>

namespace {

constexpr double kPi = 3.14159265358979323846;

// CATS3D Gaussian source convention:
//
// S(r*) = 1 / [(4*pi)^(3/2) Rout Rside Rlong]
//         * exp[-rout^2/(4 Rout^2)
//               -rside^2/(4 Rside^2)
//               -rlong^2/(4 Rlong^2)]
//
// Here R is the pair-relative Gaussian radius used in the CATS source.
double SourceGauss3D(const double rout,
                     const double rside,
                     const double rlong,
                     const double Rout,
                     const double Rside,
                     const double Rlong)
{
    if(Rout <= 0.0 || Rside <= 0.0 || Rlong <= 0.0) return 0.0;
    const double norm = 1.0 / (std::pow(4.0*kPi, 1.5)*Rout*Rside*Rlong);
    const double expo = -rout*rout/(4.0*Rout*Rout)
                        -rside*rside/(4.0*Rside*Rside)
                        -rlong*rlong/(4.0*Rlong*Rlong);
    return norm*std::exp(expo);
}

TH1D* MakeSlice(const char* name,
                const char* title,
                const char* axis,
                const unsigned nBins,
                const double rMin,
                const double rMax,
                const double Rout,
                const double Rside,
                const double Rlong)
{
    auto* h = new TH1D(name, title, nBins, rMin, rMax);
    h->SetDirectory(nullptr);

    for(unsigned i=1; i<=nBins; ++i){
        const double r = h->GetXaxis()->GetBinCenter(i);
        double rout = 0.0;
        double rside = 0.0;
        double rlong = 0.0;

        // Slice means fixing the two non-target coordinates to zero,
        // not integrating/projection over them.
        TString ax(axis);
        if(ax == "out") rout = r;
        else if(ax == "side") rside = r;
        else if(ax == "long") rlong = r;

        h->SetBinContent(i, SourceGauss3D(rout, rside, rlong,
                                          Rout, Rside, Rlong));
        h->SetBinError(i, 0.0);
    }
    return h;
}

void StyleSlice(TH1D* h, const int color)
{
    h->SetLineColor(color);
    h->SetLineWidth(3);
    h->SetMarkerSize(0);
}

void DrawOne(TH1D* h,
             const TString& outName,
             const TString& legendText,
             const bool logY)
{
    TCanvas c("cSourceSlice", "source slice", 900, 700);
    c.SetTicks(1,1);
    if(logY) c.SetLogy();

    h->SetTitle("");
    h->GetXaxis()->SetTitle("r* (fm)");
    h->GetYaxis()->SetTitle("S(r*) slice (fm^{-3})");
    h->GetXaxis()->SetTitleSize(0.052);
    h->GetYaxis()->SetTitleSize(0.052);
    h->GetXaxis()->SetLabelSize(0.045);
    h->GetYaxis()->SetLabelSize(0.045);
    h->GetYaxis()->SetTitleOffset(1.15);

    if(logY){
        h->SetMinimum(std::max(1e-8, h->GetMaximum()*1e-6));
        h->SetMaximum(h->GetMaximum()*8.0);
    }
    else{
        h->SetMinimum(0.0);
        h->SetMaximum(h->GetMaximum()*1.25);
    }

    h->Draw("HIST");

    TLegend leg(0.45, 0.76, 0.88, 0.88);
    leg.SetBorderSize(0);
    leg.SetFillStyle(0);
    leg.SetTextSize(0.04);
    leg.AddEntry(h, legendText, "l");
    leg.Draw();

    c.SaveAs(outName);
}

} // namespace

void draw_source_slice_Rfit_101(const double Rout = 4.140,
                                const double Rside = 3.837,
                                const double Rlong = 3.733,
                                const unsigned nBins = 101,
                                const double rMin = -20.0,
                                const double rMax = 20.0,
                                const char* outPrefix = "OutputFiles/source_slice_Rfit_101")
{
    gStyle->SetOptStat(0);

    std::unique_ptr<TH1D> hOut(MakeSlice("hSourceSlice_out",
                                         ";r_{out}* (fm);S(r*) slice (fm^{-3})",
                                         "out", nBins, rMin, rMax,
                                         Rout, Rside, Rlong));
    std::unique_ptr<TH1D> hSide(MakeSlice("hSourceSlice_side",
                                          ";r_{side}* (fm);S(r*) slice (fm^{-3})",
                                          "side", nBins, rMin, rMax,
                                          Rout, Rside, Rlong));
    std::unique_ptr<TH1D> hLong(MakeSlice("hSourceSlice_long",
                                          ";r_{long}* (fm);S(r*) slice (fm^{-3})",
                                          "long", nBins, rMin, rMax,
                                          Rout, Rside, Rlong));

    StyleSlice(hOut.get(), kMagenta+1);
    StyleSlice(hSide.get(), kBlue+1);
    StyleSlice(hLong.get(), kGreen+2);

    const TString prefix(outPrefix);
    TFile fout(prefix + ".root", "RECREATE");
    hOut->Write();
    hSide->Write();
    hLong->Write();
    fout.Close();

    const TString legText = Form("R=(%.3f, %.3f, %.3f) fm, %u bins in [%.0f, %.0f] fm",
                                 Rout, Rside, Rlong, nBins, rMin, rMax);

    DrawOne(hOut.get(), prefix + "_out.png", legText, false);
    DrawOne(hSide.get(), prefix + "_side.png", legText, false);
    DrawOne(hLong.get(), prefix + "_long.png", legText, false);

    // A combined plot helps compare the three source widths directly.
    TCanvas c("cSourceSliceAll", "source slices", 900, 700);
    c.SetTicks(1,1);
    hOut->SetTitle("");
    hOut->GetXaxis()->SetTitle("r* (fm)");
    hOut->GetYaxis()->SetTitle("S(r*) slice (fm^{-3})");
    hOut->GetYaxis()->SetTitleOffset(1.15);
    const double maxVal = std::max({hOut->GetMaximum(),
                                    hSide->GetMaximum(),
                                    hLong->GetMaximum()});
    hOut->SetMinimum(0.0);
    hOut->SetMaximum(maxVal*1.30);
    hOut->Draw("HIST");
    hSide->Draw("HIST SAME");
    hLong->Draw("HIST SAME");

    TLegend leg(0.52, 0.68, 0.88, 0.88);
    leg.SetBorderSize(0);
    leg.SetFillStyle(0);
    leg.SetTextSize(0.04);
    leg.AddEntry(hOut.get(), Form("out, R=%.3f fm", Rout), "l");
    leg.AddEntry(hSide.get(), Form("side, R=%.3f fm", Rside), "l");
    leg.AddEntry(hLong.get(), Form("long, R=%.3f fm", Rlong), "l");
    leg.Draw();
    c.SaveAs(prefix + "_all.png");

    std::cout << "Wrote source-slice histograms and plots with R=("
              << Rout << ", " << Rside << ", " << Rlong << ") fm\n"
              << "  " << prefix << ".root\n"
              << "  " << prefix << "_out.png\n"
              << "  " << prefix << "_side.png\n"
              << "  " << prefix << "_long.png\n"
              << "  " << prefix << "_all.png\n";
}
