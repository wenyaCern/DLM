#include "TFile.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TSystem.h"
#include <algorithm>
#include <stdexcept>

void draw_merged_nominal_slices(const char* input,
 const char* output="figures/merged_FT0C_cent1_mtMerged0_qn0_phi1",
 const char* label="FT0C plane; 10-30%; m_{T}=0.4-0.7 GeV; q_{2} bin 0, #varphi bin 1"){
 TFile f(input,"READ");if(f.IsZombie())throw std::runtime_error("Cannot open fit file");
 gStyle->SetOptStat(0);
 TCanvas c("nominal_slices","Nominal slices",1650,600);c.Divide(3,1,0.002,0.002);
 const char* dirs[]={"out","side","long"};
 for(int i=0;i<3;++i){
  c.cd(i+1);gPad->SetLeftMargin(.14);gPad->SetBottomMargin(.15);gPad->SetTopMargin(.14);
  auto* data=dynamic_cast<TH1*>(f.Get(Form("hDataCk3D_proj_%s",dirs[i])));
  auto* fit=dynamic_cast<TH1*>(f.Get(Form("hCATS3D_best_proj_%s",dirs[i])));
  auto* bkg=dynamic_cast<TH1*>(f.Get(Form("hCATS3D_background_proj_%s",dirs[i])));
  if(!data||!fit||!bkg)throw std::runtime_error("Missing unweighted projection");
  double lo=1.,hi=1.;
  for(int b=1;b<=data->GetNbinsX();++b){lo=std::min({lo,data->GetBinContent(b)-data->GetBinError(b),bkg->GetBinContent(b),fit->GetBinContent(b)});hi=std::max({hi,data->GetBinContent(b)+data->GetBinError(b),fit->GetBinContent(b)});}
  double span=std::max(.1,hi-lo);data->SetMinimum(lo-.08*span);data->SetMaximum(hi+.28*span);
  data->SetTitle(Form(";q_{%s} (GeV/c);C(q)",dirs[i]));data->SetMarkerStyle(20);data->SetMarkerSize(.7);data->SetLineColor(kBlack);
  data->GetXaxis()->SetTitleSize(.05);data->GetYaxis()->SetTitleSize(.05);data->GetYaxis()->SetTitleOffset(1.25);
  data->Draw("E1");fit->SetLineColor(kBlue);fit->SetLineWidth(2);fit->Draw("HIST SAME");
  bkg->SetLineColor(kOrange+7);bkg->SetLineStyle(2);bkg->SetLineWidth(2);bkg->Draw("HIST SAME");data->Draw("E1 SAME");
  auto* l=new TLegend(.47,.69,.93,.85);l->SetBorderSize(0);l->SetFillStyle(0);l->SetTextSize(.033);l->AddEntry(data,"Data (stored stat. errors)","pe");l->AddEntry(fit,"CATS3D fit","l");l->AddEntry(bkg,"Baseline","l");l->Draw();
  TLatex t;t.SetNDC();t.SetTextFont(42);t.SetTextSize(.030);t.DrawLatex(.14,.93,i==0 ? label : "Unweighted slice; no bootstrap band");
 }
 gSystem->mkdir(gSystem->DirName(output),true);
 c.SaveAs(Form("%s.png",output));c.SaveAs(Form("%s.pdf",output));
}
