"""Plot mixed-event weighted axial slices saved by the CATS3D fitter."""
import argparse
import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)

parser = argparse.ArgumentParser()
parser.add_argument("fit_root")
parser.add_argument("output_png")
parser.add_argument("--label", default="")
parser.add_argument("--zoom", action="store_true")
args = parser.parse_args()

source = ROOT.TFile.Open(args.fit_root)
if not source or source.IsZombie():
    raise RuntimeError(f"Cannot open {args.fit_root}")

canvas = ROOT.TCanvas("cDataVsCATS3DSlices", "", 1500, 520)
canvas.Divide(3, 1, 0.006, 0.003)
keep = []
for panel, axis in enumerate(("out", "side", "long"), 1):
    data = source.Get(f"hDataCk3D_mixedWeighted_proj_{axis}")
    fit = source.Get(f"hCATS3D_best_mixedWeighted_proj_{axis}")
    baseline = source.Get(f"hCATS3D_background_mixedWeighted_proj_{axis}")
    if not data or not fit or not baseline:
        raise RuntimeError(f"Missing {axis} projection")
    keep.extend((data, fit, baseline))
    canvas.cd(panel)
    ROOT.gPad.SetLeftMargin(0.15)
    ROOT.gPad.SetRightMargin(0.035)
    ROOT.gPad.SetBottomMargin(0.15)
    ROOT.gPad.SetTopMargin(0.14)
    ymax = max(1.27, data.GetMaximum(), fit.GetMaximum()) + 0.04
    frame = ROOT.gPad.DrawFrame(-0.10 if args.zoom else -0.20,
                               0.93 if args.zoom else 0.90,
                               0.10 if args.zoom else 0.20,
                               ymax)
    frame.SetTitle(f";q_{{{axis}}} (GeV/c);C(q)")
    frame.GetXaxis().SetTitleSize(0.052)
    frame.GetYaxis().SetTitleSize(0.052)
    frame.GetXaxis().SetLabelSize(0.043)
    frame.GetYaxis().SetLabelSize(0.043)
    frame.GetYaxis().SetTitleOffset(1.23)
    data.SetMarkerStyle(20)
    data.SetMarkerSize(0.7)
    data.SetMarkerColor(ROOT.kBlack)
    data.SetLineColor(ROOT.kBlack)
    fit.SetLineColor(ROOT.kRed + 1)
    fit.SetLineWidth(3)
    baseline.SetLineColor(ROOT.kBlue + 1)
    baseline.SetLineStyle(7)
    baseline.SetLineWidth(2)
    baseline.Draw("HIST SAME")
    fit.Draw("HIST SAME")
    data.Draw("E1 SAME")
    legend = ROOT.TLegend(0.62, 0.61, 0.96, 0.78)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.SetTextSize(0.027)
    legend.AddEntry(data, "Data", "lep")
    legend.AddEntry(fit, "CATS3D fit", "l")
    legend.AddEntry(baseline, "outer-box baseline", "l")
    legend.Draw()
    label = ROOT.TLatex()
    label.SetNDC()
    label.SetTextSize(0.035)
    label.DrawLatex(0.16, 0.92, args.label)
    label.SetTextSize(0.030)
    label.DrawLatex(0.16, 0.87, "other two |q| < 0.02 GeV/c")
    keep.extend((legend, label))

canvas.SaveAs(args.output_png)
