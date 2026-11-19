#include "TCanvas.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TLine.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TTree.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

namespace {

struct FitEntry {
    double r[3] = {0., 0., 0.};
    int status = -999;
};

struct BinSummary {
    double full[3] = {0., 0., 0.};
    double bootMean[3] = {0., 0., 0.};
    double bootSigma[3] = {0., 0., 0.};
    int nUsed = 0;
    int nMissing = 0;
    int nBadStatus = 0;
    int nBoundary = 0;
};

bool ReadFitEntry(const TString& fileName, FitEntry& out)
{
    TFile file(fileName, "READ");
    if(file.IsZombie()) return false;
    auto* tree = dynamic_cast<TTree*>(file.Get("fit_summary"));
    if(!tree || tree->GetEntries() < 1) return false;

    float rout = 0.f, rside = 0.f, rlong = 0.f, fitStatus = -999.f;
    tree->SetBranchAddress("Rout", &rout);
    tree->SetBranchAddress("Rside", &rside);
    tree->SetBranchAddress("Rlong", &rlong);
    tree->SetBranchAddress("FitStatus", &fitStatus);
    tree->GetEntry(0);
    out.r[0] = rout;
    out.r[1] = rside;
    out.r[2] = rlong;
    out.status = int(std::lround(fitStatus));
    return true;
}

double SampleSigma(const std::vector<double>& values, const double mean)
{
    if(values.size() < 2) return 0.0;
    double sum = 0.0;
    for(const double value : values) sum += (value-mean)*(value-mean);
    return std::sqrt(sum/double(values.size()-1));
}

void StyleGraph(TGraphErrors* graph, const int color, const int marker)
{
    graph->SetMarkerColor(color);
    graph->SetLineColor(color);
    graph->SetMarkerStyle(marker);
    graph->SetMarkerSize(1.25);
    graph->SetLineWidth(2);
}

} // namespace

void draw_R_vs_mT_bootstrap(const TString baseDir="farm_output")
{
    constexpr int nCent = 3;
    constexpr int nMt = 6;
    constexpr int nR = 3;
    const double mt[nMt] = {0.45, 0.55, 0.65, 0.80, 1.00, 1.305};
    const double mtErr[nMt] = {0.05, 0.05, 0.05, 0.10, 0.10, 0.205};
    const char* centLabel[nCent] = {"0#minus10%", "10#minus30%", "30#minus50%"};
    const char* rLabel[nR] = {"R_{out}", "R_{side}", "R_{long}"};
    const int colors[nR] = {kRed+1, kBlue+1, kGreen+2};
    const int markers[nR] = {20, 21, 33};

    BinSummary summary[nCent][nMt];
    bool fatal = false;
    for(int cent=0; cent<nCent; ++cent){
        for(int imt=0; imt<nMt; ++imt){
            const int mtIndex = imt + 2;
            FitEntry nominal;
            const TString nominalFile = Form(
                "%s/full/fit_CATS3D_cent%d_mt%d_nominal.root",
                baseDir.Data(), cent, mtIndex);
            if(!ReadFitEntry(nominalFile, nominal)){
                std::cerr << "Missing/invalid nominal file: " << nominalFile << "\n";
                fatal = true;
                continue;
            }
            for(int ir=0; ir<nR; ++ir) summary[cent][imt].full[ir] = nominal.r[ir];

            std::vector<double> values[nR];
            for(int replica=0; replica<100; ++replica){
                const TString bootFile = Form(
                    "%s/bootstrap/cent%d/mt%d/fit_%03d.root",
                    baseDir.Data(), cent, mtIndex, replica);
                FitEntry entry;
                if(!ReadFitEntry(bootFile, entry)){
                    ++summary[cent][imt].nMissing;
                    continue;
                }
                bool finite = true;
                bool boundary = false;
                for(int ir=0; ir<nR; ++ir){
                    finite = finite && std::isfinite(entry.r[ir]);
                    boundary = boundary || entry.r[ir] <= 1.001 || entry.r[ir] >= 15.999;
                }
                if(entry.status != 0 || !finite){
                    ++summary[cent][imt].nBadStatus;
                    continue;
                }
                if(boundary) ++summary[cent][imt].nBoundary;
                for(int ir=0; ir<nR; ++ir) values[ir].push_back(entry.r[ir]);
            }

            summary[cent][imt].nUsed = int(values[0].size());
            for(int ir=0; ir<nR; ++ir){
                double mean = 0.0;
                for(const double value : values[ir]) mean += value;
                if(!values[ir].empty()) mean /= double(values[ir].size());
                summary[cent][imt].bootMean[ir] = mean;
                summary[cent][imt].bootSigma[ir] = SampleSigma(values[ir], mean);
            }
        }
    }
    if(fatal) return;

    const TString outDir = baseDir + "/summary";
    gSystem->mkdir(outDir, true);
    std::ofstream csv((outDir + "/R_vs_mT_bootstrap.csv").Data());
    csv << "cent_index,centrality,mt_index,mt_center,mt_half_width,n_used,n_missing,n_bad_status,n_boundary,"
           "Rout_full,Rout_boot_mean,Rout_boot_sigma,Rside_full,Rside_boot_mean,Rside_boot_sigma,"
           "Rlong_full,Rlong_boot_mean,Rlong_boot_sigma\n";
    csv << std::setprecision(8);
    for(int cent=0; cent<nCent; ++cent){
        for(int imt=0; imt<nMt; ++imt){
            const auto& s = summary[cent][imt];
            csv << cent << ',' << (cent==0 ? "0-10" : cent==1 ? "10-30" : "30-50") << ','
                << imt+2 << ',' << mt[imt] << ',' << mtErr[imt] << ','
                << s.nUsed << ',' << s.nMissing << ',' << s.nBadStatus << ',' << s.nBoundary;
            for(int ir=0; ir<nR; ++ir){
                csv << ',' << s.full[ir] << ',' << s.bootMean[ir] << ',' << s.bootSigma[ir];
            }
            csv << '\n';
        }
    }
    csv.close();

    gStyle->SetOptStat(0);
    gStyle->SetEndErrorSize(5);
    gStyle->SetTitleFont(42, "XYZ");
    gStyle->SetLabelFont(42, "XYZ");
    gStyle->SetLegendFont(42);

    auto* canvas = new TCanvas("c_R_vs_mT", "R vs mT", 1650, 520);
    canvas->Divide(3, 1, 0.002, 0.002);
    TGraphErrors* graphs[nCent][nR] = {};
    for(int cent=0; cent<nCent; ++cent){
        canvas->cd(cent+1);
        gPad->SetTicks(1,1);
        gPad->SetLeftMargin(cent==0 ? 0.15 : 0.12);
        gPad->SetRightMargin(0.035);
        gPad->SetBottomMargin(0.15);
        gPad->SetTopMargin(0.09);

        double ymin = 1e9, ymax = -1e9;
        for(int imt=0; imt<nMt; ++imt) for(int ir=0; ir<nR; ++ir){
            ymin = std::min(ymin, summary[cent][imt].full[ir]-summary[cent][imt].bootSigma[ir]);
            ymax = std::max(ymax, summary[cent][imt].full[ir]+summary[cent][imt].bootSigma[ir]);
        }
        const double pad = std::max(0.25, 0.18*(ymax-ymin));

        for(int ir=0; ir<nR; ++ir){
            double y[nMt], ey[nMt];
            for(int imt=0; imt<nMt; ++imt){
                y[imt] = summary[cent][imt].full[ir];
                ey[imt] = summary[cent][imt].bootSigma[ir];
            }
            graphs[cent][ir] = new TGraphErrors(nMt, mt, y, mtErr, ey);
            graphs[cent][ir]->SetName(Form("g_%s_cent%d", ir==0?"Rout":ir==1?"Rside":"Rlong", cent));
            StyleGraph(graphs[cent][ir], colors[ir], markers[ir]);
            if(ir==0){
                graphs[cent][ir]->SetTitle(Form("%s;#it{m}_{T} (GeV/#it{c}^{2});R (fm)", centLabel[cent]));
                graphs[cent][ir]->GetXaxis()->SetLimits(0.34, 1.53);
                graphs[cent][ir]->GetYaxis()->SetRangeUser(std::max(0.0,ymin-pad), ymax+pad);
                graphs[cent][ir]->GetXaxis()->SetTitleSize(0.055);
                graphs[cent][ir]->GetYaxis()->SetTitleSize(0.055);
                graphs[cent][ir]->GetXaxis()->SetLabelSize(0.045);
                graphs[cent][ir]->GetYaxis()->SetLabelSize(0.045);
                graphs[cent][ir]->Draw("AP");
            } else {
                graphs[cent][ir]->Draw("P SAME");
            }
        }
        auto* legend = new TLegend(0.57, 0.64, 0.93, 0.88);
        legend->SetBorderSize(0);
        legend->SetFillStyle(0);
        legend->SetTextSize(0.045);
        for(int ir=0; ir<nR; ++ir) legend->AddEntry(graphs[cent][ir], rLabel[ir], "pe");
        legend->Draw();
        auto* label = new TLatex(0.18, 0.85, Form("Centrality %s", centLabel[cent]));
        label->SetNDC();
        label->SetTextFont(42);
        label->SetTextSize(0.050);
        label->Draw();
    }

    canvas->SaveAs(outDir + "/R_vs_mT_bootstrap.png");
    canvas->SaveAs(outDir + "/R_vs_mT_bootstrap.pdf");

    const int centColors[nCent] = {kRed+1, kBlue+1, kGreen+2};
    const int centMarkers[nCent] = {20, 21, 33};
    TCanvas* componentCanvas[nR] = {};
    TGraphErrors* componentGraphs[nR][nCent] = {};
    for(int ir=0; ir<nR; ++ir){
        const char* componentName = ir==0 ? "Rout" : ir==1 ? "Rside" : "Rlong";
        componentCanvas[ir] = new TCanvas(Form("c_%s_vs_mT", componentName),
                                          Form("%s vs mT", componentName), 760, 650);
        componentCanvas[ir]->SetTicks(1,1);
        componentCanvas[ir]->SetLeftMargin(0.14);
        componentCanvas[ir]->SetRightMargin(0.04);
        componentCanvas[ir]->SetBottomMargin(0.13);
        componentCanvas[ir]->SetTopMargin(0.06);

        double ymin = 1e9, ymax = -1e9;
        for(int cent=0; cent<nCent; ++cent) for(int imt=0; imt<nMt; ++imt){
            ymin = std::min(ymin, summary[cent][imt].full[ir]-summary[cent][imt].bootSigma[ir]);
            ymax = std::max(ymax, summary[cent][imt].full[ir]+summary[cent][imt].bootSigma[ir]);
        }
        const double pad = std::max(0.25, 0.12*(ymax-ymin));

        for(int cent=0; cent<nCent; ++cent){
            double y[nMt], ey[nMt];
            for(int imt=0; imt<nMt; ++imt){
                y[imt] = summary[cent][imt].full[ir];
                ey[imt] = summary[cent][imt].bootSigma[ir];
            }
            componentGraphs[ir][cent] = new TGraphErrors(nMt, mt, y, mtErr, ey);
            componentGraphs[ir][cent]->SetName(Form("g_%s_cent%d_comparison", componentName, cent));
            StyleGraph(componentGraphs[ir][cent], centColors[cent], centMarkers[cent]);
            if(cent==0){
                componentGraphs[ir][cent]->SetTitle(Form(";#it{m}_{T} (GeV/#it{c}^{2});%s (fm)", rLabel[ir]));
                componentGraphs[ir][cent]->GetXaxis()->SetLimits(0.34, 1.53);
                componentGraphs[ir][cent]->GetYaxis()->SetRangeUser(std::max(0.0,ymin-pad), ymax+pad);
                componentGraphs[ir][cent]->GetXaxis()->SetTitleSize(0.052);
                componentGraphs[ir][cent]->GetYaxis()->SetTitleSize(0.052);
                componentGraphs[ir][cent]->GetXaxis()->SetLabelSize(0.043);
                componentGraphs[ir][cent]->GetYaxis()->SetLabelSize(0.043);
                componentGraphs[ir][cent]->Draw("AP");
            } else {
                componentGraphs[ir][cent]->Draw("P SAME");
            }
        }
        auto* legend = new TLegend(0.61, 0.68, 0.92, 0.91);
        legend->SetBorderSize(0);
        legend->SetFillStyle(0);
        legend->SetTextSize(0.040);
        for(int cent=0; cent<nCent; ++cent){
            legend->AddEntry(componentGraphs[ir][cent], Form("%s", centLabel[cent]), "pe");
        }
        legend->Draw();
        componentCanvas[ir]->SaveAs(outDir + Form("/%s_vs_mT_centrality.png", componentName));
        componentCanvas[ir]->SaveAs(outDir + Form("/%s_vs_mT_centrality.pdf", componentName));
    }

    TFile rootOut(outDir + "/R_vs_mT_bootstrap.root", "RECREATE");
    canvas->Write();
    for(int cent=0; cent<nCent; ++cent) for(int ir=0; ir<nR; ++ir) graphs[cent][ir]->Write();
    for(int ir=0; ir<nR; ++ir){
        componentCanvas[ir]->Write();
        for(int cent=0; cent<nCent; ++cent) componentGraphs[ir][cent]->Write();
    }
    rootOut.Close();

    std::cout << "Wrote " << outDir << "/R_vs_mT_bootstrap.{png,pdf,root}\n";
    std::cout << "Wrote " << outDir << "/R_vs_mT_bootstrap.csv\n";
    for(int cent=0; cent<nCent; ++cent){
        for(int imt=0; imt<nMt; ++imt){
            const auto& s = summary[cent][imt];
            std::cout << "cent" << cent << " mt" << imt+2
                      << " used=" << s.nUsed << " missing=" << s.nMissing
                      << " bad=" << s.nBadStatus << " boundary=" << s.nBoundary
                      << " Rout=" << s.full[0] << "+/-" << s.bootSigma[0]
                      << " Rside=" << s.full[1] << "+/-" << s.bootSigma[1]
                      << " Rlong=" << s.full[2] << "+/-" << s.bootSigma[2] << '\n';
        }
    }
}
