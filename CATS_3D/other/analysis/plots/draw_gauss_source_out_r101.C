#include "TCanvas.h"
#include "TFile.h"
#include "TH3D.h"
#include "TH1D.h"
#include "TLatex.h"
#include "TParameter.h"
#include "TSystem.h"
#include "TStyle.h"
#include <stdexcept>
R__LOAD_LIBRARY(build/libFitCATS3D)
double Basics_Source_Gauss3D(double*);

void draw_gauss_source_out_r101(double Rout=5.,double Rside=5.,double Rlong=5.){
 if(Rout<=0||Rside<=0||Rlong<=0)throw std::runtime_error("Radii must be positive");
 gStyle->SetOptStat(0);
 const int n=101;const double lo=-20.,hi=20.;
 TH3D source("hSourceGauss3D",";r_{out}^{*} (fm);r_{side}^{*} (fm);r_{long}^{*} (fm)",n,lo,hi,n,lo,hi,n,lo,hi);
 source.SetDirectory(nullptr);
 double pars[9]={0,0,0,0,0,0,Rout,Rside,Rlong};
 for(int i=1;i<=n;++i)for(int j=1;j<=n;++j)for(int k=1;k<=n;++k){
  pars[3]=source.GetXaxis()->GetBinCenter(i);pars[4]=source.GetYaxis()->GetBinCenter(j);pars[5]=source.GetZaxis()->GetBinCenter(k);
  source.SetBinContent(i,j,k,Basics_Source_Gauss3D(pars));
 }
 TH1D slice("hSourceOutCentralSlice",";r_{out}^{*} (fm);S(r_{out}^{*},0,0) (fm^{-3})",n,lo,hi);
 slice.SetDirectory(nullptr);
 for(int i=1;i<=n;++i){slice.SetBinContent(i,source.GetBinContent(i,51,51));slice.SetBinError(i,0);}
 TCanvas c("cSourceOut","",1050,750);c.SetLeftMargin(.18);c.SetBottomMargin(.16);c.SetTopMargin(.07);c.SetRightMargin(.04);c.SetTicks(1,1);
 slice.GetXaxis()->SetTitleSize(.043);slice.GetYaxis()->SetTitleSize(.043);
 slice.GetXaxis()->SetLabelSize(.038);slice.GetYaxis()->SetLabelSize(.038);
 slice.GetXaxis()->SetTitleOffset(1.15);
 slice.SetLineColor(kViolet+1);slice.SetLineWidth(3);slice.SetMinimum(0);slice.SetMaximum(slice.GetMaximum()*1.32);slice.GetYaxis()->SetTitleOffset(1.65);slice.Draw("HIST");
 TLatex lab;lab.SetNDC();lab.SetTextFont(42);lab.SetTextSize(.034);
 lab.DrawLatex(.18,.87,Form("Gaussian source: R_{out/side/long} = (%.2f, %.2f, %.2f) fm",Rout,Rside,Rlong));
 lab.DrawLatex(.18,.82,"Central slice: r_{side}^{*} = r_{long}^{*} = 0");
 lab.DrawLatex(.18,.77,"101^{3} grid, each axis [-20, 20] fm; #Deltar^{*} = 0.39604 fm");
 gSystem->mkdir("figures",true);
 TString stem=Form("figures/gauss_source_out_r101_R%.2f_%.2f_%.2f",Rout,Rside,Rlong);
 c.SaveAs(stem+".png");c.SaveAs(stem+".pdf");
 TFile f(stem+".root","RECREATE");source.Write();slice.Write();c.Write();
 TParameter<double>("Rout_fm",Rout).Write();TParameter<double>("Rside_fm",Rside).Write();TParameter<double>("Rlong_fm",Rlong).Write();
 printf("Grid step %.9f fm; central slice peak %.10g fm^-3; source volume integral %.9f\n",40./101,slice.GetBinContent(51),source.Integral("width"));
}
