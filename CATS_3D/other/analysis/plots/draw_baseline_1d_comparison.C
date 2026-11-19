#include "TCanvas.h"
#include "TH1.h"
#include "TF1.h"
#include "TBox.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TArrow.h"

// Exaggerated conceptual curves, not fits to data. The continuous-region
// example assumes far-region information dominates; not a general property.
void draw_baseline_1d_comparison(){
 gStyle->SetOptStat(0);gSystem->mkdir("figures",true);
 TCanvas combined("baseline_email","Baseline selection analogy",2400,760);
 combined.Divide(2,1,0.002,0.002);
 for(int comparison=0;comparison<2;++comparison){
  double ymin=comparison?.83:.94;
  TCanvas c(Form("analogy%d",comparison),"1D illustration",1200,760);
  c.SetLeftMargin(.10);c.SetRightMargin(.035);c.SetBottomMargin(.13);c.SetTopMargin(.07);
  auto* frame=c.DrawFrame(0.,ymin,1.,1.63);frame->SetTitle(";q;C(q)");
  for(auto* ax:{frame->GetXaxis(),frame->GetYaxis()}){ax->SetLabelSize(0);ax->SetTickLength(0);ax->SetTitleSize(.052);}
  frame->GetXaxis()->SetTitleOffset(.60);frame->GetYaxis()->SetTitleOffset(.75);
  TBox near(.26,ymin,.43,1.63),far(.76,ymin,.96,1.63);
  for(auto* b:{&near,&far}){b->SetFillColorAlpha(kViolet+1,.15);b->SetLineColor(0);b->Draw();}
  TF1 cf(Form("cf%d",comparison),"(1.02+0.07*x*x)*(1+0.46*exp(-pow(x/0.095,2)))+0.007*sin(145*x+0.3)+0.005*sin(237*x+1.1)+0.003*sin(391*x)",0,1);
  cf.SetLineColor(kBlack);cf.SetLineWidth(3);cf.SetNpx(1500);cf.Draw("SAME");
  TF1 base("base","1.02+0.07*x*x",0,1),bn("bn","1.02+0.07*x*x",.26,.43),bf("bf","1.02+0.07*x*x",.76,.96);
  base.SetLineColor(kBlue+1);base.SetLineStyle(7);base.SetLineWidth(3);
  bn.SetLineColor(kBlue+1);bf.SetLineColor(kBlue+1);bn.SetLineWidth(8);bf.SetLineWidth(8);
  TF1 red("near_only","1.02+0.07*x*x-1.2*(x-0.26)*(x-0.43)",0,.43);
  TF1 green("far_only","1.02+0.07*x*x+0.12*pow((0.85-x)/0.55,2)",0,1);
  TF1 magenta("continuous","1.02+0.07*x*x+0.44*(x-0.62)*(x-0.92)",0,1);
  TLegend leg(.44,comparison?.60:.65,.76,comparison?.84:.78);
  leg.SetBorderSize(0);leg.SetFillStyle(0);leg.SetTextSize(.032);leg.AddEntry(&cf,"CF","l");
  if(!comparison){base.Draw("SAME");bn.Draw("SAME");bf.Draw("SAME");leg.AddEntry(&bn,"Near + far baseline","l");}
  else{
   red.SetNpx(1000);red.SetLineColor(kRed+1);red.SetLineWidth(4);red.Draw("SAME");
   green.SetLineColor(kGreen+2);green.SetLineWidth(4);green.Draw("SAME");
   magenta.SetLineColor(kMagenta);magenta.SetLineWidth(4);magenta.SetLineStyle(7);magenta.Draw("SAME");
   leg.AddEntry(&red,"Near-only baseline","l");leg.AddEntry(&green,"Far-only baseline","l");leg.AddEntry(&magenta,"Continuous-region baseline","l");
  }
  leg.Draw();TLatex t;t.SetTextFont(42);t.SetTextAlign(22);t.SetTextColor(kViolet+1);t.SetTextSize(.033);
  t.DrawLatex(.345,1.57,"Near sideband");t.DrawLatex(.86,1.57,"Far sideband");
  t.SetNDC();t.SetTextAlign(12);t.SetTextColor(kGray+2);t.SetTextSize(.023);
  t.DrawLatex(.10,.965,"Schematic");
  TArrow deviation(.34,1.035,.34,1.115,.015,"<>");
  if(comparison){deviation.SetLineColor(kGreen+2);deviation.SetFillColor(kGreen+2);deviation.SetLineWidth(2);deviation.Draw();}
  c.SaveAs(comparison?"figures/baseline_1d_sideband_alternatives.png":"figures/baseline_1d_two_sidebands_minimal.png");
  combined.cd(comparison+1);c.DrawClonePad();
 }
 combined.SaveAs("figures/baseline_1d_email_comparison.png");
}
