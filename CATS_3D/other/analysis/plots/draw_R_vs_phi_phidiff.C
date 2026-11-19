#include "TCanvas.h"
#include "TFile.h"
#include "TGraphErrors.h"
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

struct FitEntryPhi {
    double r[3] = {0.,0.,0.};
    int status = -999;
};

bool ReadFitPhi(const TString& fileName, FitEntryPhi& out)
{
    TFile file(fileName, "READ");
    if(file.IsZombie()) return false;
    auto* tree = dynamic_cast<TTree*>(file.Get("fit_summary"));
    if(!tree || tree->GetEntries()<1) return false;
    float rout=0.f, rside=0.f, rlong=0.f, status=-999.f;
    tree->SetBranchAddress("Rout", &rout);
    tree->SetBranchAddress("Rside", &rside);
    tree->SetBranchAddress("Rlong", &rlong);
    tree->SetBranchAddress("FitStatus", &status);
    tree->GetEntry(0);
    out.r[0]=rout; out.r[1]=rside; out.r[2]=rlong;
    out.status=int(std::lround(status));
    return true;
}

double SampleSigmaPhi(const std::vector<double>& values)
{
    if(values.size()<2) return 0.0;
    double mean=0.0;
    for(double value: values) mean+=value;
    mean/=double(values.size());
    double sum=0.0;
    for(double value: values) sum+=(value-mean)*(value-mean);
    return std::sqrt(sum/double(values.size()-1));
}

} // namespace

void draw_R_vs_phi_phidiff(const TString baseDir="farm_output")
{
    constexpr int nCent=3, nMt=6, nQn=3, nPhi=6, nR=3;
    // Exact merged-bin geometry from the input THnSparse phi axis:
    // 12 bins over [-0.2,3.4] rad, merged pairwise.
    const double phi[nPhi] = {0.1,0.7,1.3,1.9,2.5,3.1};
    const double phiErr[nPhi] = {0.3,0.3,0.3,0.3,0.3,0.3};
    const double mtLow[nMt]  = {0.40,0.50,0.60,0.70,0.90,1.10};
    const double mtHigh[nMt] = {0.50,0.60,0.70,0.90,1.10,1.51};
    const char* centText[nCent] = {"0#minus10%","10#minus30%","30#minus50%"};
    const char* centCsv[nCent] = {"0-10","10-30","30-50"};
    const char* rName[nR] = {"Rout","Rside","Rlong"};
    const char* rTitle[nR] = {"R_{out}","R_{side}","R_{long}"};
    const int qnColor[nQn] = {kRed+1,kBlue+1,kGreen+2};
    const int qnMarker[nQn] = {20,21,33};
    const double splitScale = std::sqrt(18.0);

    double fullSigma[nCent][nMt][nR] = {};
    int fullUsed[nCent][nMt] = {};
    bool fatal=false;
    for(int cent=0;cent<nCent;++cent) for(int imt=0;imt<nMt;++imt){
        const int mt=imt+2;
        std::vector<double> values[nR];
        for(int replica=0;replica<100;++replica){
            FitEntryPhi entry;
            const TString fn=Form("%s/bootstrap/cent%d/mt%d/fit_%03d.root",
                                  baseDir.Data(),cent,mt,replica);
            if(!ReadFitPhi(fn,entry) || entry.status!=0) continue;
            bool finite=true;
            for(int ir=0;ir<nR;++ir) finite=finite&&std::isfinite(entry.r[ir]);
            if(!finite) continue;
            for(int ir=0;ir<nR;++ir) values[ir].push_back(entry.r[ir]);
        }
        fullUsed[cent][imt]=int(values[0].size());
        if(fullUsed[cent][imt]<2){
            std::cerr<<"Insufficient bootstrap fits for cent"<<cent<<" mt"<<mt<<"\n";
            fatal=true;
        }
        for(int ir=0;ir<nR;++ir) fullSigma[cent][imt][ir]=SampleSigmaPhi(values[ir]);
    }
    if(fatal) return;

    FitEntryPhi result[nCent][nMt][nQn][nPhi];
    int missing=0,badStatus=0,boundary=0;
    for(int cent=0;cent<nCent;++cent) for(int imt=0;imt<nMt;++imt)
      for(int qn=0;qn<nQn;++qn) for(int iphi=0;iphi<nPhi;++iphi){
        const int mt=imt+2;
        const TString fn=Form("%s/phidiff/cent%d/mt%d/qn%d/fit_phi%d.root",
                              baseDir.Data(),cent,mt,qn,iphi);
        if(!ReadFitPhi(fn,result[cent][imt][qn][iphi])){
            ++missing; std::cerr<<"Missing/invalid phi-diff file: "<<fn<<"\n"; continue;
        }
        if(result[cent][imt][qn][iphi].status!=0) ++badStatus;
        for(int ir=0;ir<nR;++ir){
            const double value=result[cent][imt][qn][iphi].r[ir];
            if(value<=1.001||value>=15.999) {++boundary; break;}
        }
      }
    if(missing){
        std::cerr<<"Abort plotting: missing/invalid phi-diff files="<<missing<<"\n";
        return;
    }

    const TString outDir=baseDir+"/phidiff_summary";
    gSystem->mkdir(outDir,true);
    std::ofstream csv((outDir+"/R_vs_phi_phidiff.csv").Data());
    csv<<"cent_index,centrality,mt_index,mt_low,mt_high,qn_index,phi_index,phi_low,phi_high,phi_center,"
          "Rout,Rout_sigma_est,Rside,Rside_sigma_est,Rlong,Rlong_sigma_est,fit_status\n";
    csv<<std::setprecision(8);
    for(int cent=0;cent<nCent;++cent) for(int imt=0;imt<nMt;++imt)
      for(int qn=0;qn<nQn;++qn) for(int iphi=0;iphi<nPhi;++iphi){
        const auto& e=result[cent][imt][qn][iphi];
        csv<<cent<<','<<centCsv[cent]<<','<<imt+2<<','<<mtLow[imt]<<','<<mtHigh[imt]<<','
           <<qn<<','<<iphi<<','<<phi[iphi]-phiErr[iphi]<<','<<phi[iphi]+phiErr[iphi]<<','<<phi[iphi];
        for(int ir=0;ir<nR;++ir) csv<<','<<e.r[ir]<<','<<splitScale*fullSigma[cent][imt][ir];
        csv<<','<<e.status<<'\n';
      }
    csv.close();

    gStyle->SetOptStat(0);
    gStyle->SetEndErrorSize(5);
    TFile rootOut(outDir+"/R_vs_phi_phidiff.root","RECREATE");
    for(int cent=0;cent<nCent;++cent){
      const TString centDir=outDir+Form("/cent%d",cent);
      gSystem->mkdir(centDir,true);
      for(int imt=0;imt<nMt;++imt){
        const int mt=imt+2;
        auto* canvas=new TCanvas(Form("c_R_vs_phi_cent%d_mt%d",cent,mt),"",1650,520);
        canvas->Divide(3,1,0.002,0.002);
        for(int ir=0;ir<nR;++ir){
          canvas->cd(ir+1);
          gPad->SetTicks(1,1);
          gPad->SetLeftMargin(ir==0?0.15:0.12);
          gPad->SetRightMargin(0.035);
          gPad->SetBottomMargin(0.16);
          gPad->SetTopMargin(0.08);
          double ymin=1e9,ymax=-1e9;
          const double ey=splitScale*fullSigma[cent][imt][ir];
          for(int qn=0;qn<nQn;++qn) for(int iphi=0;iphi<nPhi;++iphi){
            const double y=result[cent][imt][qn][iphi].r[ir];
            ymin=std::min(ymin,y-ey); ymax=std::max(ymax,y+ey);
          }
          const double pad=std::max(0.20,0.16*(ymax-ymin));
          TGraphErrors* graphs[nQn]={};
          for(int qn=0;qn<nQn;++qn){
            double y[nPhi],yerr[nPhi];
            for(int iphi=0;iphi<nPhi;++iphi){
              y[iphi]=result[cent][imt][qn][iphi].r[ir];
              yerr[iphi]=ey;
            }
            graphs[qn]=new TGraphErrors(nPhi,phi,y,phiErr,yerr);
            graphs[qn]->SetName(Form("g_%s_cent%d_mt%d_qn%d",rName[ir],cent,mt,qn));
            graphs[qn]->SetMarkerStyle(qnMarker[qn]);
            graphs[qn]->SetMarkerColor(qnColor[qn]);
            graphs[qn]->SetLineColor(qnColor[qn]);
            graphs[qn]->SetMarkerSize(1.2);
            graphs[qn]->SetLineWidth(2);
            if(qn==0){
              graphs[qn]->SetTitle(Form(";%s;%s (fm)","#varphi_{pair} #minus #Psi_{EP} (rad)",rTitle[ir]));
              graphs[qn]->GetXaxis()->SetLimits(-0.22,3.42);
              graphs[qn]->GetYaxis()->SetRangeUser(std::max(0.0,ymin-pad),ymax+pad);
              graphs[qn]->GetXaxis()->SetTitleSize(0.050);
              graphs[qn]->GetYaxis()->SetTitleSize(0.055);
              graphs[qn]->GetXaxis()->SetLabelSize(0.043);
              graphs[qn]->GetYaxis()->SetLabelSize(0.043);
              graphs[qn]->Draw("AP");
            } else graphs[qn]->Draw("P SAME");
            rootOut.cd(); graphs[qn]->Write();
          }
          auto* legend=new TLegend(0.68,0.66,0.93,0.89);
          legend->SetBorderSize(0); legend->SetFillStyle(0); legend->SetTextSize(0.043);
          for(int qn=0;qn<nQn;++qn) legend->AddEntry(graphs[qn],Form("q_{n} bin %d",qn),"pe");
          legend->Draw();
          auto* label=new TLatex(0.18,0.86,Form("%s, %.2f < #it{m}_{T} < %.2f",centText[cent],mtLow[imt],mtHigh[imt]));
          label->SetNDC(); label->SetTextFont(42); label->SetTextSize(0.045); label->Draw();
        }
        const TString stem=centDir+Form("/R_vs_phi_cent%d_mt%d",cent,mt);
        canvas->SaveAs(stem+".png");
        canvas->SaveAs(stem+".pdf");
        rootOut.cd(); canvas->Write();
      }
    }
    rootOut.Close();
    std::cout<<"phi-diff quality: files=324 missing="<<missing
             <<" badStatus="<<badStatus<<" boundary="<<boundary<<"\n";
    std::cout<<"uncertainty prescription: sigma(phi,qn)=sqrt(18)*sigma(full bootstrap)\n";
    std::cout<<"wrote "<<outDir<<"\n";
}
