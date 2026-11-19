#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TStyle.h"

void draw_corrected_plane_closure(){
 const double x[3]={0,1,2},ex[3]={0,0,0};
 const double eCobs[3]={.0370693,.0651506,.0865815},rC[3]={.545802134,.812512719,.908091740};
 const double eAobs[3]={.0684452,.0887440,.1031710},rA[3]={.603096070,.666918642,.734383966};
 const double sCobs[3]={.0256570,.0256638,.0255409},sAobs[3]={.0258137,.0258333,.0257231};
 double yC[3],yA[3],eC[3],eA[3];for(int i=0;i<3;++i){yC[i]=eCobs[i]/rC[i];yA[i]=eAobs[i]/rA[i];eC[i]=sCobs[i]/rC[i];eA[i]=sAobs[i]/rA[i];}
 gStyle->SetOptStat(0);auto*c=new TCanvas("c_plane_closure","",760,650);c->SetLeftMargin(.14);c->SetBottomMargin(.13);c->SetTicks();
 auto*gC=new TGraphErrors(3,x,yC,ex,eC);auto*gA=new TGraphErrors(3,x,yA,ex,eA);
 gC->SetTitle(";q_{2}^{FT0C} bin;#epsilon_{final}/Res_{2}(#Psi_{2})");gC->SetMarkerStyle(20);gC->SetMarkerColor(kBlue+1);gC->SetLineColor(kBlue+1);gC->SetLineWidth(3);gC->SetMarkerSize(1.3);gC->GetYaxis()->SetRangeUser(0.,.20);gC->Draw("APL");
 gA->SetMarkerStyle(21);gA->SetMarkerColor(kRed+1);gA->SetLineColor(kRed+1);gA->SetLineWidth(3);gA->SetMarkerSize(1.3);gA->Draw("PL SAME");
 auto*l=new TLegend(.19,.72,.57,.89);l->SetBorderSize(0);l->SetFillStyle(0);l->AddEntry(gC,"#Psi_{2}^{FT0C}, corrected by Res(FT0C)","pe");l->AddEntry(gA,"#Psi_{2}^{FT0A}, corrected by Res(FT0A)","pe");l->Draw();
 c->SaveAs("farm_output/plane_resolution_closure_cent1_mt0p4to0p7.png");c->SaveAs("farm_output/plane_resolution_closure_cent1_mt0p4to0p7.pdf");
 for(int i=0;i<3;++i)printf("q%d C=%.6f +/- %.6f A=%.6f +/- %.6f\n",i,yC[i],eC[i],yA[i],eA[i]);
}
