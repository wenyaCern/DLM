#include "TCanvas.h"
#include "TH1.h"
#include "TF1.h"
#include "TBox.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TSystem.h"

// Conceptual illustration only: not measured data or a literal 3D range mapping.
void draw_baseline_1d_analogy(){
 gStyle->SetOptStat(0);
 TCanvas c("baseline_analogy","1D baseline sideband illustration",1200,780);
 c.SetLeftMargin(.11);c.SetRightMargin(.035);c.SetBottomMargin(.22);c.SetTopMargin(.13);
 auto* frame=c.DrawFrame(0.,.96,1.,1.65);
 frame->SetTitle(";Relative momentum q (schematic);C(q)");
 frame->GetXaxis()->SetLabelSize(0);frame->GetXaxis()->SetTickLength(0);
 frame->GetXaxis()->SetTitleSize(.045);frame->GetXaxis()->SetTitleOffset(.9);
 frame->GetYaxis()->SetTitleSize(.048);frame->GetYaxis()->SetTitleOffset(1.0);
 frame->GetYaxis()->SetLabelSize(0);frame->GetYaxis()->SetTickLength(0);
 TBox femto(0.,.96,.23,1.65);femto.SetFillColorAlpha(kGray,.14);femto.SetLineColor(0);femto.Draw();
 TBox near(.26,.96,.43,1.65);near.SetFillColorAlpha(kViolet+1,.17);near.SetLineColor(0);near.Draw();
 TBox far(.76,.96,.96,1.65);far.SetFillColorAlpha(kOrange+1,.23);far.SetLineColor(0);far.Draw();
 TF1 cf("illustrative_cf","(1.02+0.07*x*x)*(1+0.46*exp(-pow(x/0.095,2))) + 0.007*sin(145*x+0.3) + 0.005*sin(237*x+1.1) + 0.003*sin(391*x)",0,1);
 cf.SetLineColor(kBlack);cf.SetLineWidth(3);cf.SetNpx(1000);cf.Draw("SAME");
 TF1 baseline("illustrative_baseline","1.02+0.07*x*x",0,1);
 baseline.SetLineColor(kBlue+1);baseline.SetLineStyle(7);baseline.SetLineWidth(3);baseline.Draw("SAME");
 TF1 bn("baseline_near","1.02+0.07*x*x",.26,.43);bn.SetLineColor(kBlue+1);bn.SetLineWidth(8);bn.Draw("SAME");
 TF1 bf("baseline_far","1.02+0.07*x*x",.76,.96);bf.SetLineColor(kBlue+1);bf.SetLineWidth(8);bf.Draw("SAME");
 TLatex t;t.SetNDC();t.SetTextFont(42);t.SetTextColor(kBlack);t.SetTextSize(.037);
 t.DrawLatex(.11,.95,"Why use two baseline sidebands?  A 1D analogy");
 t.SetTextSize(.025);t.DrawLatex(.11,.90,"Conceptual illustration; not a literal projection of the 3D hybrid selection");
 t.SetTextAlign(22);t.SetTextSize(.027);
 t.DrawLatex(.205,.82,"Femto region");t.DrawLatex(.205,.775,"Excluded from prefit");
 t.SetTextColor(kViolet+1);t.DrawLatex(.402,.82,"Near sideband");t.DrawLatex(.402,.775,"Outside the peak");
 t.SetTextColor(kOrange+7);t.DrawLatex(.84,.82,"Far sideband");t.DrawLatex(.84,.775,"Long-range constraint");
 t.SetTextAlign(12);t.SetTextColor(kBlack);
 TLegend leg(.48,.57,.73,.70);leg.SetBorderSize(0);leg.SetFillStyle(0);leg.SetTextSize(.029);
 leg.AddEntry(&cf,"CF (illustrative)","l");leg.AddEntry(&baseline,"Baseline model","l");leg.Draw();
 t.SetTextSize(.027);t.DrawLatex(.39,.38,"Dashed: baseline extrapolation / interpolation");
 t.SetTextSize(.028);t.DrawLatex(.11,.125,"Near + far sidebands jointly constrain the baseline level and curvature.");
 t.DrawLatex(.11,.075,"The baseline is then extrapolated beneath the femtoscopic peak.");
 c.RedrawAxis();gSystem->mkdir("figures",true);
 c.SaveAs("figures/baseline_1d_two_sidebands_analogy.png");
 c.SaveAs("figures/baseline_1d_two_sidebands_analogy.pdf");
}
