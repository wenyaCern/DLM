#include "Basics3D.h"
#include "/Users/wuwenya/Documents/code/femto_flow/CATS/CATS_install/install/include/CATSconstants.h"

#include "TCanvas.h"
#include "TFile.h"
#include "TH1.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TStyle.h"

#include <array>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
void Configure(CATS3D& kitty,int along,const std::array<double,3>& radius){
    for(int axis=0;axis<3;++axis){
        const unsigned n=(axis==along)?10:2;
        const double hi=(axis==along)?100.0:20.0;
        if(axis==0)kitty.SetMomBinsOut(n,0.0,hi);
        if(axis==1)kitty.SetMomBinsSide(n,0.0,hi);
        if(axis==2)kitty.SetMomBinsLong(n,0.0,hi);
    }
    kitty.SetGridX(31,-40.0,40.0);
    kitty.SetGridY(31,-40.0,40.0);
    kitty.SetGridZ(31,-40.0,40.0);
    kitty.SetFastKernelMode(true);
    kitty.SetFastKernelSubdivisions(1);
    CATSparameters pars(CATSparameters::tSource3D,3,true);
    for(unsigned a=0;a<3;++a)pars.SetParameter(a,radius[a]);
    kitty.SetAnaSource(Basics_Source_Gauss3D,pars);
    kitty.SetPdgId(211,211);
    kitty.SetQuantumStatistics(1);
    kitty.SetQ1Q2(1);
    kitty.SetRedMass(0.5*Mass_Pic);
    kitty.SetWfMode(CATS3D::wfCoulomb);
}
}

int main(int argc,char** argv){
    const std::string choice=(argc>1)?argv[1]:"085";
    if(choice!="083"&&choice!="085")throw std::runtime_error("Use 083 or 085");
    const double lambdaGen=(choice=="083")?0.83:0.85;
    const double lambda=lambdaGen*(1.0-0.14)*(1.0-0.14);
    const std::array<double,3> radius=(choice=="083")
        ?std::array<double,3>{8.98989,7.99913,9.40678}
        :std::array<double,3>{9.09009,8.08106,9.50016};
    const TString directory="derivedFull_Bminus_cent0to10_anyCPR_fineKernel_derivedLambda_baseline100to200_fit/mt0/qn0/";
    const TString fitPath=directory+Form("fit_phi0_lambda%s_rinit8p0.root",choice.c_str());
    const TString outBase=directory+Form("approx_r40_fixedR_lambda%s",choice.c_str());
    TFile fit(fitPath,"READ");
    if(fit.IsZombie())throw std::runtime_error("Cannot open fitted data");
    TFile out(outBase+".root","RECREATE");
    gStyle->SetOptStat(0);
    TCanvas canvas("cApproxR40FixedR","",1500,550);
    canvas.Divide(3,1,.006,.003);
    std::vector<std::unique_ptr<TH1>> keepApprox;
    const char* axisNames[]={"out","side","long"};
    for(int axis=0;axis<3;++axis){
        const TString tag=axisNames[axis];
        auto* data=dynamic_cast<TH1*>(fit.Get("hDataCk3D_mixedWeighted_proj_"+tag));
        auto* old=dynamic_cast<TH1*>(fit.Get("hCATS3D_best_mixedWeighted_proj_"+tag));
        auto* baseline=dynamic_cast<TH1*>(fit.Get("hCATS3D_background_mixedWeighted_proj_"+tag));
        if(!data||!old||!baseline)throw std::runtime_error("Missing fit projection");
        std::unique_ptr<TH1> approx(dynamic_cast<TH1*>(baseline->Clone("hApproxR40_"+tag)));
        approx->SetDirectory(nullptr);
        CATS3D kitty;
        Configure(kitty,axis,radius);
        kitty.KillTheCat();
        for(int ibin=0;ibin<5;++ibin){
            double avg=0.0;
            for(int dx=0;dx<2;++dx)
                for(int dy=0;dy<2;++dy)
                    for(int dz=0;dz<2;++dz){
                        unsigned i[3]={unsigned(dx),unsigned(dy),unsigned(dz)};
                        i[axis]=unsigned(2*ibin+i[axis]);
                        avg+=(1.0+lambda*(kitty.GetCorrFun(i[0],i[1],i[2])-1.0))/8.0;
                    }
            const double q=0.01+0.02*ibin;
            for(double sign:{-1.0,1.0}){
                const int bin=approx->FindBin(sign*q);
                approx->SetBinContent(bin,avg*baseline->GetBinContent(bin));
            }
            if(ibin<3)
                std::cout<<tag<<" q="<<q<<" data="<<data->GetBinContent(data->FindBin(q))
                         <<" r20fit="<<old->GetBinContent(old->FindBin(q))
                         <<" r40fixed="<<approx->GetBinContent(approx->FindBin(q))<<'\n';
        }
        out.cd();approx->Write();
        canvas.cd(axis+1);
        gPad->SetLeftMargin(.15);gPad->SetRightMargin(.035);
        gPad->SetBottomMargin(.15);gPad->SetTopMargin(.15);
        auto* frame=gPad->DrawFrame(-.10,.92,.10,1.36);
        frame->SetTitle(Form(";q_{%s} (GeV/c);C(q)",axisNames[axis]));
        frame->GetYaxis()->SetTitleOffset(1.15);
        old->SetLineColor(kRed+1);old->SetLineWidth(2);old->SetLineStyle(2);
        approx->SetLineColor(kBlue+1);approx->SetLineWidth(3);
        data->SetMarkerStyle(20);data->SetMarkerSize(.65);
        data->SetLineColor(kBlack);
        old->Draw("HIST SAME");approx->Draw("HIST SAME");data->Draw("E1 SAME");
        if(axis==0){
            auto* leg=new TLegend(.49,.68,.95,.84);
            leg->SetBorderSize(0);leg->SetFillStyle(0);leg->SetTextSize(.032);
            leg->AddEntry(data,"Data","lep");
            leg->AddEntry(old,"r = #pm20 fm fit","l");
            leg->AddEntry(approx.get(),"r = #pm40 fm, fixed R","l");
            leg->Draw();
        }
        TLatex title;title.SetNDC();title.SetTextSize(.031);
        title.DrawLatex(.16,.92,Form("0-10%%, 0.2<m_{T}<0.4; #lambda_{gen}=%.2f",lambdaGen));
        title.SetTextSize(.027);
        title.DrawLatex(.16,.87,"blue: approximate, no refit");
        keepApprox.push_back(std::move(approx));
    }
    out.cd();canvas.Write();out.Close();
    canvas.SaveAs(TString(outBase)+".png");
}
