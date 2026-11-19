#include "TFile.h"
#include "TH1.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TSystem.h"
#include <stdexcept>

void draw_merged_lambda_comparison(){
 TFile f("File/lam/lambda_merged_cent10to30.root","READ");
 const char* tags[]={"0.40to0.70","0.70to1.50"};
 const char* labels[]={"0.4 #leq m_{T} < 0.7 GeV","0.7 #leq m_{T} < 1.5 GeV"};
 TGraph* g[2];
 double low=1.,high=0.;
 for(int i=0;i<2;++i){
  auto* h=dynamic_cast<TH1*>(f.Get(Form("hPurityKstar_mT%s_cent10to30",tags[i])));
  if(!h)throw std::runtime_error("Missing merged lambda histogram");
  g[i]=new TGraph();g[i]->SetName(Form("lambdaInput_mt%d",i));
  for(int b=1;b<=h->GetNbinsX();++b){
   double x=h->GetBinCenter(b)*1000.;if(x>200.)break;
   double y=h->GetBinContent(b);
   g[i]->SetPoint(g[i]->GetN(),x,y);low=std::min(low,y);high=std::max(high,y);
  }
  g[i]->SetLineColor(i==0?kViolet+1:kOrange+7);g[i]->SetMarkerColor(i==0?kViolet+1:kOrange+7);
  g[i]->SetLineWidth(3);g[i]->SetMarkerStyle(i==0?20:21);g[i]->SetMarkerSize(.85);
 }
 gStyle->SetOptStat(0);TCanvas c("lambda_comparison","Merged lambda",900,680);
 c.SetLeftMargin(.14);c.SetBottomMargin(.14);
 auto* frame=c.DrawFrame(0.,std::max(0.,low-.04),200.,std::min(1.,high+.12));
 frame->SetTitle(";k* (MeV/c);#lambda_{input}(k*)");
 frame->GetXaxis()->SetTitleSize(.045);frame->GetYaxis()->SetTitleSize(.045);
 frame->GetYaxis()->SetTitleOffset(1.35);
 for(auto* graph:g)graph->Draw("LP SAME");
 TLegend leg(.49,.70,.88,.85);leg.SetBorderSize(0);leg.SetFillStyle(0);
 for(int i=0;i<2;++i)leg.AddEntry(g[i],labels[i],"lp");leg.Draw();
 TLatex t;t.SetNDC();t.SetTextFont(42);t.SetTextSize(.035);
 t.DrawLatex(.17,.86,"Pb-Pb, 10-30%");
 t.SetTextSize(.030);t.DrawLatex(.17,.79,"Before weak-decay scaling");
 gSystem->mkdir("figures",true);
 c.SaveAs("figures/lambda_merged_mt_input_comparison.png");c.SaveAs("figures/lambda_merged_mt_input_comparison.pdf");
 TFile out("figures/lambda_merged_mt_input_comparison.root","RECREATE");for(auto* graph:g)graph->Write();c.Write();
}
