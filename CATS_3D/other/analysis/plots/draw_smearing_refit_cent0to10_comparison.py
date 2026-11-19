"""Compare the inclusive-PID and 0-10% smearing CATS3D fits in axial slices."""

from pathlib import Path
import ROOT


ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
directory = Path(__file__).resolve().parents[3] / (
    "fit/results/derivedFull_Bminus_cent0to10_anyCPR_adaptiveKernel_smear3d_weak155_fit"
    "/mt0/qn0"
)
old = ROOT.TFile.Open(str(directory / "fit_phi0_PIDresponse_diagonalBand_fixed_kernelFitBoxOnly.root"))
new = ROOT.TFile.Open(str(directory / "fit_phi0_PIDresponse_diagonalBand_fixed_kernelFitBoxOnly_cent0to10Smear.root"))
if not old or old.IsZombie() or not new or new.IsZombie():
    raise RuntimeError("Missing old or new fit ROOT file")

canvas = ROOT.TCanvas("cSmearingRefit", "", 1500, 530)
canvas.Divide(3, 1)
keep = []
for panel, axis in enumerate(("out", "side", "long"), 1):
    canvas.cd(panel)
    ROOT.gPad.SetLeftMargin(0.15)
    ROOT.gPad.SetRightMargin(0.035)
    ROOT.gPad.SetBottomMargin(0.15)
    ROOT.gPad.SetTopMargin(0.14)
    data = new.Get(f"hDataCk3D_mixedWeighted_proj_{axis}")
    old_fit = old.Get(f"hCATS3D_best_mixedWeighted_proj_{axis}")
    new_fit = new.Get(f"hCATS3D_best_mixedWeighted_proj_{axis}")
    baseline = new.Get(f"hCATS3D_background_mixedWeighted_proj_{axis}")
    if not all((data, old_fit, new_fit, baseline)):
        raise RuntimeError(f"Missing projection for {axis}")
    frame = ROOT.gPad.DrawFrame(-0.2, 0.95, 0.2, 1.28)
    frame.SetTitle(f";q_{{{axis}}} (GeV/c);C(q)")
    frame.GetXaxis().SetTitleSize(0.05)
    frame.GetYaxis().SetTitleSize(0.05)
    frame.GetYaxis().SetTitleOffset(1.26)
    data.SetMarkerStyle(20)
    data.SetMarkerSize(0.72)
    data.SetLineColor(ROOT.kBlack)
    old_fit.SetLineColor(ROOT.kBlue + 1)
    old_fit.SetLineStyle(7)
    old_fit.SetLineWidth(3)
    new_fit.SetLineColor(ROOT.kRed + 1)
    new_fit.SetLineWidth(3)
    baseline.SetLineColor(ROOT.kGray + 2)
    baseline.SetLineStyle(3)
    baseline.SetLineWidth(2)
    baseline.Draw("HIST SAME")
    old_fit.Draw("HIST SAME")
    new_fit.Draw("HIST SAME")
    data.Draw("E1 SAME")
    legend = ROOT.TLegend(0.53, 0.58, 0.96, 0.81)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.SetTextSize(0.027)
    legend.AddEntry(data, "Data", "lep")
    legend.AddEntry(old_fit, "inclusive-cent smearing", "l")
    legend.AddEntry(new_fit, "0-10% smearing", "l")
    legend.AddEntry(baseline, "baseline", "l")
    legend.Draw()
    title = ROOT.TLatex()
    title.SetNDC()
    title.SetTextSize(0.031)
    title.DrawLatex(0.16, 0.92, "0-10%, 0.2 < m_{T} < 0.4; q_{2}/#phi bin 0")
    title.DrawLatex(0.16, 0.87, "other two |q| < 0.02 GeV/c")
    keep.extend((data, old_fit, new_fit, baseline, frame, legend, title))

output = directory / "fit_phi0_inclusive_vs_cent0to10_smearing_slices.png"
canvas.SaveAs(str(output))
print(output)
