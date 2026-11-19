#include "TCanvas.h"
#include "TH1.h"
#include "TF1.h"
#include "TBox.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TSystem.h"

// Conceptual illustration only; horizontal coordinates are not physical cuts.
void draw_normalization_1d_schematic() {
 gStyle->SetOptStat(0);
 gSystem->mkdir("figures",true);
 TCanvas c("normalization_schematic","1D normalization illustration",1200,760);
 c.SetLeftMargin(.10);c.SetRightMargin(.035);
 c.SetBottomMargin(.13);c.SetTopMargin(.07);
 auto* frame=c.DrawFrame(0.,.94,1.,1.63);
 frame->SetTitle(";q;C(q)");
 for(auto* ax:{frame->GetXaxis(),frame->GetYaxis()}) {
  ax->SetLabelSize(0);ax->SetTickLength(0);ax->SetTitleSize(.052);
 }
 frame->GetXaxis()->SetTitleOffset(.60);
 frame->GetYaxis()->SetTitleOffset(.75);
 TBox region(.58,.94,.96,1.63);
 region.SetFillColorAlpha(kViolet+1,.18);
 region.SetLineColor(0);region.Draw();
 TF1 cf("normalization_cf","(1.02+0.07*x*x)*(1+0.46*exp(-pow(x/0.095,2)))+0.007*sin(145*x+0.3)+0.005*sin(237*x+1.1)+0.003*sin(391*x)",0,1);
 cf.SetLineColor(kBlack);cf.SetLineWidth(3);cf.SetNpx(1500);cf.Draw("SAME");
 TLatex text;text.SetTextFont(42);text.SetTextAlign(22);
 text.SetTextColor(kViolet+1);text.SetTextSize(.035);
 text.DrawLatex(.77,1.56,"Normalization region");
 text.SetNDC();text.SetTextAlign(12);text.SetTextColor(kGray+2);text.SetTextSize(.023);
 text.DrawLatex(.10,.965,"Schematic");
 c.SaveAs("figures/normalization_1d_schematic.png");
}
