#include "TCanvas.h"
#include "TFile.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TTree.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

namespace {
bool readR(const TString& fn,double r[3]){
  TFile f(fn,"READ"); auto*t=dynamic_cast<TTree*>(f.Get("fit_summary"));
  if(f.IsZombie()||!t||!t->GetEntries())return false;
  float ro,rs,rl;t->SetBranchAddress("Rout",&ro);t->SetBranchAddress("Rside",&rs);
  t->SetBranchAddress("Rlong",&rl);t->GetEntry(0);r[0]=ro;r[1]=rs;r[2]=rl;return true;
}
double rms(const std::vector<double>&v){double a=0,s=0;for(double x:v)a+=x;a/=v.size();for(double x:v)s+=(x-a)*(x-a);return std::sqrt(s/(v.size()-1));}
void linfit(const double*x,const double*y,const double*e,int n,double&A,double&B,double&eA,double&eB,double&cov){
  double s=0,sx=0,sxx=0,sy=0,sxy=0;for(int i=0;i<n;++i){double z=std::cos(2*x[i]),w=1/(e[i]*e[i]);s+=w;sx+=w*z;sxx+=w*z*z;sy+=w*y[i];sxy+=w*z*y[i];}
  double d=s*sxx-sx*sx;A=(sxx*sy-sx*sxy)/d;B=(s*sxy-sx*sy)/d;eA=std::sqrt(sxx/d);eB=std::sqrt(s/d);cov=-sx/d;
}
}

void draw_ft0a_rough_q2(const TString base="farm_output",int cent=1,int mt=3,
                        bool merged047=false,const TString plane="ft0a"){
  const double phi[6]={.1,.7,1.3,1.9,2.5,3.1}; const char*cn[3]={"out","side","long"};
  const int col[3]={kRed+1,kBlue+1,kGreen+2},mk[3]={20,21,33}; double sig[3]={};
  for(int b=0;b<100;++b){double r[3];if(!readR(Form("%s/bootstrap/cent%d/mt%d/fit_%03d.root",base.Data(),cent,mt,b),r))continue;for(int j=0;j<3;++j){/* collected below */}}
  for(int j=0;j<3;++j){
    double invVar=0.;
    const int first=merged047?2:mt, last=merged047?4:mt;
    for(int im=first;im<=last;++im){std::vector<double>v;for(int b=0;b<100;++b){double r[3];if(readR(Form("%s/bootstrap/cent%d/mt%d/fit_%03d.root",base.Data(),cent,im,b),r))v.push_back(r[j]);}const double s=rms(v);if(s>0)invVar+=1/(s*s);}
    sig[j]=(invVar>0?1/std::sqrt(invVar):0)*std::sqrt(18.);
  }
  double eps[3]={},eeps[3]={}; auto*c=new TCanvas("c_ft0a_rough","",1650,520);c->Divide(3,1,.002,.002);
  for(int j=0;j<3;++j){c->cd(j+1);gPad->SetTicks();double lo=1e9,hi=-1e9;TGraphErrors*g[3];TGraph*f[3];
    for(int q=0;q<3;++q){double y[6],ey[6],ex[6];for(int p=0;p<6;++p){double r[3];const TString fn=merged047?Form("%s/phidiff_%s_merged047/qn%d/fit_phi%d.root",base.Data(),plane.Data(),q,p):Form("%s/phidiff_ft0a/cent%d/mt%d/qn%d/fit_phi%d.root",base.Data(),cent,mt,q,p);readR(fn,r);y[p]=r[j]*r[j];ey[p]=2*r[j]*sig[j];ex[p]=.3;lo=std::min(lo,y[p]-ey[p]);hi=std::max(hi,y[p]+ey[p]);}
      double A,B,eA,eB,cov;linfit(phi,y,ey,6,A,B,eA,eB,cov);if(j==1){eps[q]=B/A;eeps[q]=std::sqrt(eB*eB/(A*A)+B*B*eA*eA/std::pow(A,4)-2*B*cov/std::pow(A,3));}
      g[q]=new TGraphErrors(6,phi,y,ex,ey);g[q]->SetMarkerStyle(mk[q]);g[q]->SetMarkerColor(col[q]);g[q]->SetLineColor(col[q]);
      double xx[201],yy[201];for(int i=0;i<201;++i){xx[i]=-.2+3.6*i/200.;yy[i]=A+B*std::cos(2*xx[i]);}f[q]=new TGraph(201,xx,yy);f[q]->SetLineColor(col[q]);f[q]->SetLineWidth(3);
      if(q==0){g[q]->SetTitle(Form(";#Delta#varphi=#varphi_{pair}-#Psi_{2}^{%s} (rad);R_{%s}^{2} (fm^{2})",plane=="ft0c"?"FT0C":"FT0A",cn[j]));g[q]->GetXaxis()->SetLimits(-.2,3.4);g[q]->GetYaxis()->SetRangeUser(lo-.1*(hi-lo),hi+.1*(hi-lo));g[q]->Draw("AP");}else g[q]->Draw("P SAME");f[q]->Draw("L SAME");}
    auto*l=new TLegend(.66,.68,.92,.9);l->SetBorderSize(0);for(int q=0;q<3;++q)l->AddEntry(g[q],Form("q_{2} bin %d",q),"pe");l->Draw();
  }
  TString od=base+(merged047?Form("/phidiff_%s_summary/rough_cent1_mt0p4to0p7",plane.Data()):"/phidiff_ft0a_summary/rough_cent1_mt3");gSystem->mkdir(od,true);c->SaveAs(od+"/R2_vs_phi.png");c->SaveAs(od+"/R2_vs_phi.pdf");
  auto*ce=new TCanvas("c_eps_ft0a_rough","",700,600);double x[3]={0,1,2},ex[3]={0,0,0};auto*ge=new TGraphErrors(3,x,eps,ex,eeps);ge->SetTitle(";q_{2} bin;#epsilon_{final}=2R_{side,2}^{2}/R_{side,0}^{2}");ge->SetMarkerStyle(20);ge->SetMarkerSize(1.4);ge->SetLineWidth(2);ge->Draw("APL");ce->SaveAs(od+"/epsilon_vs_q2.png");ce->SaveAs(od+"/epsilon_vs_q2.pdf");
  for(int q=0;q<3;++q)std::cout<<"q2 "<<q<<" epsilon="<<eps[q]<<" +/- "<<eeps[q]<<"\n";
}
