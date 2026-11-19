#include "TCanvas.h"
#include "TFile.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TTree.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

namespace {
struct EFitEntry { double r[3]={0.,0.,0.}; int status=-999; };
struct FourierResult {
  double A=0., B=0., eA=0., eB=0., covAB=0., chi2=0.;
  int ndf=0;
  double Epsilon() const { return A!=0. ? B/A : 0.; }
  double EpsilonError() const {
    if(A==0.) return 0.;
    const double variance=eB*eB/(A*A)+B*B*eA*eA/(A*A*A*A)-2.*B*covAB/(A*A*A);
    return std::sqrt(std::max(0.,variance));
  }
};

bool ReadE(const TString& fn, EFitEntry& out){
  TFile f(fn,"READ"); if(f.IsZombie()) return false;
  auto*t=dynamic_cast<TTree*>(f.Get("fit_summary")); if(!t||t->GetEntries()<1)return false;
  float ro=0,rs=0,rl=0,st=-999;
  t->SetBranchAddress("Rout",&ro); t->SetBranchAddress("Rside",&rs);
  t->SetBranchAddress("Rlong",&rl); t->SetBranchAddress("FitStatus",&st); t->GetEntry(0);
  out.r[0]=ro;out.r[1]=rs;out.r[2]=rl;out.status=int(std::lround(st)); return true;
}
double SampleSigmaE(const std::vector<double>& v){
  if(v.size()<2)return 0.; double m=0.;for(double x:v)m+=x;m/=v.size();
  double s=0.;for(double x:v)s+=(x-m)*(x-m);return std::sqrt(s/double(v.size()-1));
}
FourierResult FitCos2(const double* phi,const double* y,const double* ey,const int n){
  // Weighted linear least squares in the basis {1, cos(2 phi)}.
  double s=0.,sx=0.,sxx=0.,sy=0.,sxy=0.;
  for(int i=0;i<n;++i){
    const double x=std::cos(2.*phi[i]); const double w=1./(ey[i]*ey[i]);
    s+=w;sx+=w*x;sxx+=w*x*x;sy+=w*y[i];sxy+=w*x*y[i];
  }
  const double det=s*sxx-sx*sx;
  FourierResult r;
  r.A=(sxx*sy-sx*sxy)/det; r.B=(s*sxy-sx*sy)/det;
  r.eA=std::sqrt(sxx/det); r.eB=std::sqrt(s/det); r.covAB=-sx/det;
  for(int i=0;i<n;++i){
    const double d=(y[i]-r.A-r.B*std::cos(2.*phi[i]))/ey[i]; r.chi2+=d*d;
  }
  r.ndf=n-2; return r;
}
}

void draw_epsilon_vs_q2(const TString baseDir="farm_output",
                        const TString phiFitDir="phidiff",
                        const TString summaryDir="phidiff_summary"){
  constexpr int NC=3,NM=6,NQ=3,NP=6,NR=3;
  const double phi[NP]={0.1,0.7,1.3,1.9,2.5,3.1};
  const double mtLo[NM]={.40,.50,.60,.70,.90,1.10};
  const double mtHi[NM]={.50,.60,.70,.90,1.10,1.51};
  const double mtCenter[NM]={.45,.55,.65,.80,1.00,1.305};
  const double mtHalfWidth[NM]={.05,.05,.05,.10,.10,.205};
  const char* centLabel[NC]={"0#minus10%","10#minus30%","30#minus50%"};
  const char* centCsv[NC]={"0-10","10-30","30-50"};
  const char* comp[NR]={"out","side","long"};
  const int colors[NM]={kRed+1,kOrange+7,kGreen+2,kCyan+2,kBlue+1,kViolet+1};
  const int markers[NM]={20,21,22,23,33,34};
  const int qnColor[NQ]={kRed+1,kBlue+1,kGreen+2};
  const int qnMarker[NQ]={20,21,33};
  const double scale=std::sqrt(18.0);

  double sigmaFull[NC][NM][NR]={};
  for(int c=0;c<NC;++c)for(int m=0;m<NM;++m){
    std::vector<double> v[NR];
    for(int b=0;b<100;++b){ EFitEntry e;
      if(!ReadE(Form("%s/bootstrap/cent%d/mt%d/fit_%03d.root",baseDir.Data(),c,m+2,b),e)||e.status!=0)continue;
      for(int r=0;r<NR;++r)v[r].push_back(e.r[r]);
    }
    for(int r=0;r<NR;++r)sigmaFull[c][m][r]=SampleSigmaE(v[r]);
  }

  FourierResult fit[NC][NM][NQ][NR];
  EFitEntry point[NC][NM][NQ][NP];
  int missing=0,bad=0;
  for(int c=0;c<NC;++c)for(int m=0;m<NM;++m)for(int q=0;q<NQ;++q){
    for(int p=0;p<NP;++p){
      if(!ReadE(Form("%s/%s/cent%d/mt%d/qn%d/fit_phi%d.root",baseDir.Data(),phiFitDir.Data(),c,m+2,q,p),point[c][m][q][p]))++missing;
      else if(point[c][m][q][p].status!=0)++bad;
    }
    for(int r=0;r<NR;++r){
      double y[NP],ey[NP];
      for(int p=0;p<NP;++p){
        const double R=point[c][m][q][p].r[r];
        y[p]=R*R; ey[p]=2.*std::abs(R)*scale*sigmaFull[c][m][r];
      }
      fit[c][m][q][r]=FitCos2(phi,y,ey,NP);
    }
  }
  if(missing){std::cerr<<"Missing phi-diff files: "<<missing<<"\n";return;}

  const TString outDir=baseDir+"/"+summaryDir;gSystem->mkdir(outDir,true);
  std::ofstream csv((outDir+"/R2_fourier_epsilon_vs_q2.csv").Data());
  csv<<"cent_index,centrality,mt_index,mt_low,mt_high,q2_bin,component,R2_0,R2_0_err,cos_coefficient_B,B_err,R2_2,cov_A_B,chi2,ndf,epsilon_side,epsilon_side_err\n";
  csv<<std::setprecision(9);
  for(int c=0;c<NC;++c)for(int m=0;m<NM;++m)for(int q=0;q<NQ;++q)for(int r=0;r<NR;++r){
    const auto& f=fit[c][m][q][r];
    csv<<c<<','<<centCsv[c]<<','<<m+2<<','<<mtLo[m]<<','<<mtHi[m]<<','<<q<<','<<comp[r]<<','
       <<f.A<<','<<f.eA<<','<<f.B<<','<<f.eB<<','<<f.B/2.<<','<<f.covAB<<','<<f.chi2<<','<<f.ndf<<',';
    if(r==1)csv<<f.Epsilon()<<','<<f.EpsilonError(); else csv<<',';
    csv<<'\n';
  } csv.close();

  gStyle->SetOptStat(0);gStyle->SetEndErrorSize(5);
  TFile fout(outDir+"/epsilon_vs_q2.root","RECREATE");

  const TString r2Dir=outDir+"/R2_vs_phi";
  gSystem->mkdir(r2Dir,true);
  for(int c=0;c<NC;++c){
    const TString centDir=r2Dir+Form("/cent%d",c);
    gSystem->mkdir(centDir,true);
    for(int m=0;m<NM;++m){
      auto* cr2=new TCanvas(Form("c_R2_vs_phi_cent%d_mt%d",c,m+2),"",1650,520);
      cr2->Divide(3,1,.002,.002);
      for(int r=0;r<NR;++r){
        cr2->cd(r+1);gPad->SetTicks(1,1);gPad->SetBottomMargin(.16);
        gPad->SetLeftMargin(r==0?.15:.12);gPad->SetRightMargin(.025);gPad->SetTopMargin(.08);
        double ymin=1e99,ymax=-1e99;
        for(int q=0;q<NQ;++q)for(int p=0;p<NP;++p){
          const double R=point[c][m][q][p].r[r];
          const double y=R*R, ey=2.*std::abs(R)*scale*sigmaFull[c][m][r];
          ymin=std::min(ymin,y-ey);ymax=std::max(ymax,y+ey);
        }
        const double pad=std::max(1.0,.14*(ymax-ymin));
        TGraphErrors* gd[NQ]={}; TGraph* gf[NQ]={};
        for(int q=0;q<NQ;++q){
          double y[NP],ey[NP],ex[NP]={.3,.3,.3,.3,.3,.3};
          for(int p=0;p<NP;++p){const double R=point[c][m][q][p].r[r];y[p]=R*R;ey[p]=2.*std::abs(R)*scale*sigmaFull[c][m][r];}
          gd[q]=new TGraphErrors(NP,phi,y,ex,ey);
          gd[q]->SetName(Form("g_R2_%s_cent%d_mt%d_q%d",comp[r],c,m+2,q));
          gd[q]->SetMarkerStyle(qnMarker[q]);gd[q]->SetMarkerColor(qnColor[q]);gd[q]->SetLineColor(qnColor[q]);gd[q]->SetMarkerSize(1.15);gd[q]->SetLineWidth(1);
          if(q==0){gd[q]->SetTitle(Form(";#varphi_{pair} #minus #Psi_{EP} (rad);R_{%s}^{2} (fm^{2})",comp[r]));gd[q]->GetXaxis()->SetLimits(-.22,3.42);gd[q]->GetYaxis()->SetRangeUser(std::max(0.,ymin-pad),ymax+pad);gd[q]->GetXaxis()->SetTitleSize(.050);gd[q]->GetYaxis()->SetTitleSize(.052);gd[q]->GetXaxis()->SetLabelSize(.043);gd[q]->GetYaxis()->SetLabelSize(.043);gd[q]->Draw("AP");}else gd[q]->Draw("P SAME");
          constexpr int NCURVE=241; double xc[NCURVE],yc[NCURVE];
          for(int i=0;i<NCURVE;++i){xc[i]=-.2+3.6*double(i)/double(NCURVE-1);yc[i]=fit[c][m][q][r].A+fit[c][m][q][r].B*std::cos(2.*xc[i]);}
          gf[q]=new TGraph(NCURVE,xc,yc);gf[q]->SetName(Form("gFit_R2_%s_cent%d_mt%d_q%d",comp[r],c,m+2,q));gf[q]->SetLineColor(qnColor[q]);gf[q]->SetLineWidth(3);gf[q]->Draw("L SAME");
          fout.cd();gd[q]->Write();gf[q]->Write();
        }
        auto* leg=new TLegend(.65,.64,.94,.89);leg->SetBorderSize(0);leg->SetFillStyle(0);leg->SetTextSize(.040);
        for(int q=0;q<NQ;++q)leg->AddEntry(gd[q],Form("q_{2} bin %d",q),"pe");leg->Draw();
        auto* lab=new TLatex(.18,.86,Form("%s, %.2f < #it{m}_{T} < %.2f",centLabel[c],mtLo[m],mtHi[m]));lab->SetNDC();lab->SetTextFont(42);lab->SetTextSize(.045);lab->Draw();
      }
      const TString stem=centDir+Form("/R2_vs_phi_cent%d_mt%d",c,m+2);
      cr2->SaveAs(stem+".png");cr2->SaveAs(stem+".pdf");fout.cd();cr2->Write();
    }
  }

  auto* canvas=new TCanvas("c_epsilon_vs_q2","",1650,520);canvas->Divide(3,1,.002,.002);
  for(int c=0;c<NC;++c){canvas->cd(c+1);gPad->SetTicks(1,1);gPad->SetBottomMargin(.15);gPad->SetLeftMargin(.17);gPad->SetRightMargin(.025);gPad->SetTopMargin(.08);
    double ymin=1e9,ymax=-1e9;
    for(int m=0;m<NM;++m)for(int q=0;q<NQ;++q){const auto& f=fit[c][m][q][1];ymin=std::min(ymin,f.Epsilon()-f.EpsilonError());ymax=std::max(ymax,f.Epsilon()+f.EpsilonError());}
    const double pad=std::max(.02,.15*(ymax-ymin));
    auto* frame=new TH1D(Form("hframe_epsilon_cent%d",c),"",3,-.5,2.5);
    frame->SetDirectory(nullptr);
    frame->SetMinimum(ymin-pad); frame->SetMaximum(ymax+pad);
    frame->GetXaxis()->SetTitle("q_{2} bin");
    frame->GetYaxis()->SetTitle("#epsilon_{final}=2R_{side,2}^{2}/R_{side,0}^{2}");
    frame->GetXaxis()->SetBinLabel(1,"0"); frame->GetXaxis()->SetBinLabel(2,"1"); frame->GetXaxis()->SetBinLabel(3,"2");
    frame->GetXaxis()->SetTitleSize(.052); frame->GetYaxis()->SetTitleSize(.050);
    frame->GetXaxis()->SetLabelSize(.048); frame->GetYaxis()->SetLabelSize(.043);
    frame->GetYaxis()->SetTitleOffset(1.45);
    frame->Draw("AXIS");
    TGraphErrors* gr[NM]={};
    for(int m=0;m<NM;++m){double x[NQ],y[NQ],ex[NQ]={0.,0.,0.},ey[NQ];
      for(int q=0;q<NQ;++q){x[q]=q+0.035*(m-2.5);const auto& f=fit[c][m][q][1];y[q]=f.Epsilon();ey[q]=f.EpsilonError();}
      gr[m]=new TGraphErrors(NQ,x,y,ex,ey);gr[m]->SetName(Form("g_epsilon_cent%d_mt%d",c,m+2));gr[m]->SetMarkerColor(colors[m]);gr[m]->SetLineColor(colors[m]);gr[m]->SetMarkerStyle(markers[m]);gr[m]->SetMarkerSize(1.15);gr[m]->SetLineWidth(2);
      gr[m]->Draw("PL SAME");fout.cd();gr[m]->Write();
    }
    auto* leg=new TLegend(.58,.58,.94,.90);leg->SetBorderSize(0);leg->SetFillStyle(0);leg->SetTextSize(.035);
    for(int m=0;m<NM;++m)leg->AddEntry(gr[m],Form("%.2f#minus%.2f",mtLo[m],mtHi[m]),"pe");leg->Draw();
    auto* lab=new TLatex(.18,.87,Form("Centrality %s",centLabel[c]));lab->SetNDC();lab->SetTextFont(42);lab->SetTextSize(.047);lab->Draw();
  }
  canvas->SaveAs(outDir+"/epsilon_vs_q2.png");canvas->SaveAs(outDir+"/epsilon_vs_q2.pdf");fout.cd();canvas->Write();

  auto* cMt=new TCanvas("c_epsilon_vs_mT","",1650,520);cMt->Divide(3,1,.002,.002);
  for(int c=0;c<NC;++c){
    cMt->cd(c+1);gPad->SetTicks(1,1);gPad->SetBottomMargin(.15);gPad->SetLeftMargin(.17);gPad->SetRightMargin(.025);gPad->SetTopMargin(.08);
    double ymin=1e9,ymax=-1e9;
    for(int m=0;m<NM;++m)for(int q=0;q<NQ;++q){const auto& f=fit[c][m][q][1];ymin=std::min(ymin,f.Epsilon()-f.EpsilonError());ymax=std::max(ymax,f.Epsilon()+f.EpsilonError());}
    const double pad=std::max(.02,.15*(ymax-ymin));
    TGraphErrors* gr[NQ]={};
    for(int q=0;q<NQ;++q){double y[NM],ey[NM];
      for(int m=0;m<NM;++m){const auto& f=fit[c][m][q][1];y[m]=f.Epsilon();ey[m]=f.EpsilonError();}
      gr[q]=new TGraphErrors(NM,mtCenter,y,mtHalfWidth,ey);
      gr[q]->SetName(Form("g_epsilon_vs_mT_cent%d_q%d",c,q));gr[q]->SetMarkerColor(qnColor[q]);gr[q]->SetLineColor(qnColor[q]);gr[q]->SetMarkerStyle(qnMarker[q]);gr[q]->SetMarkerSize(1.15);gr[q]->SetLineWidth(2);
      if(q==0){gr[q]->SetTitle(";#it{m}_{T} (GeV/#it{c}^{2});#epsilon_{final}=2R_{side,2}^{2}/R_{side,0}^{2}");gr[q]->GetXaxis()->SetLimits(.34,1.53);gr[q]->GetYaxis()->SetRangeUser(ymin-pad,ymax+pad);gr[q]->GetXaxis()->SetTitleSize(.052);gr[q]->GetYaxis()->SetTitleSize(.050);gr[q]->GetXaxis()->SetLabelSize(.043);gr[q]->GetYaxis()->SetLabelSize(.043);gr[q]->GetYaxis()->SetTitleOffset(1.45);gr[q]->Draw("APL");}else gr[q]->Draw("PL SAME");
      fout.cd();gr[q]->Write();
    }
    auto* leg=new TLegend(.66,.69,.94,.90);leg->SetBorderSize(0);leg->SetFillStyle(0);leg->SetTextSize(.041);
    for(int q=0;q<NQ;++q)leg->AddEntry(gr[q],Form("q_{2} bin %d",q),"pe");leg->Draw();
    auto* lab=new TLatex(.18,.87,Form("Centrality %s",centLabel[c]));lab->SetNDC();lab->SetTextFont(42);lab->SetTextSize(.047);lab->Draw();
  }
  cMt->SaveAs(outDir+"/epsilon_vs_mT.png");cMt->SaveAs(outDir+"/epsilon_vs_mT.pdf");fout.cd();cMt->Write();fout.Close();
  std::cout<<"quality: phi fits=324 missing="<<missing<<" badStatus="<<bad<<"\n";
  std::cout<<"wrote "<<outDir<<"/epsilon_vs_q2.{png,pdf,root}\n";
  std::cout<<"wrote "<<outDir<<"/epsilon_vs_mT.{png,pdf}\n";
  std::cout<<"wrote "<<outDir<<"/R2_fourier_epsilon_vs_q2.csv\n";
}
