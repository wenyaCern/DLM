"""Compare fixed-R toy response before/after a 5% decrease in Rside."""
import csv
import math
import sys

import numpy as np
import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
base = "fit/results/derivedFull_Bminus_cent0to10_anyCPR_fineKernel_derivedLambda_baseline100to200_fit/mt0/qn0/"
fit = ROOT.TFile.Open(base + "fit_phi0_lambda083_rinit8p0.root")
all_mode = len(sys.argv) > 1 and sys.argv[1] == "all"
scale = (0.84 / 0.86) ** 2  # omega 0.14 -> 0.16
mapping = np.array([9-i if i < 10 else i-10 for i in range(20)])


def load_model(filename):
    octant = np.empty((10, 10, 10), dtype=float)
    with open(base + filename, newline="") as stream:
        for row in csv.DictReader(stream):
            octant[int(row["i_out"]), int(row["i_side"]), int(row["i_long"])] = float(row["C"])
    return 1 + (octant[np.ix_(mapping, mapping, mapping)] - 1) * scale


def blur_axis(array, sigma_mev, axis):
    sigma = sigma_mev / 5.0
    radius = math.ceil(4 * sigma)
    offsets = np.arange(-radius, radius + 1)
    weights = np.exp(-0.5 * (offsets / sigma) ** 2)
    weights /= weights.sum()
    padding = [(0, 0)] * 3
    padding[axis] = (radius, radius)
    padded = np.pad(array, padding, mode="constant", constant_values=1.0)
    output = np.zeros_like(array)
    for offset, weight in zip(offsets, weights):
        slices = [slice(None)] * 3
        slices[axis] = slice(radius + offset, radius + offset + array.shape[axis])
        output += weight * padded[tuple(slices)]
    return output


def smear(array):
    for axis, sigma in enumerate((8.0, 2.6, 2.6)):
        array = blur_axis(array, sigma, axis)
    return array


original = smear(load_model("cats3d_r40_lambda083_q5_octant.csv"))
changed = smear(load_model(
    "cats3d_r40_lambda083_q5_Rall95_octant.csv" if all_mode
    else "cats3d_r40_lambda083_q5_Rside95_octant.csv"))
models = {"original R": original,
          "all R x 0.95" if all_mode else "Rside x 0.95": changed}
if all_mode:
    peak_raw = changed[10:14, 6:14, 6:14].mean()
    hdata = fit.Get("hDataCk3D_mixedWeighted_proj_out")
    hbase = fit.Get("hCATS3D_background_mixedWeighted_proj_out")
    qbin = hdata.FindBin(0.01)
    target_raw = hdata.GetBinContent(qbin) / hbase.GetBinContent(hbase.FindBin(0.01))
    lambda_scale = (target_raw - 1.0) / (peak_raw - 1.0)
    print("peak-aligned lambda_gen", 0.83 * lambda_scale)
    models["all R x 0.95, peak aligned"] = 1.0 + (changed - 1.0) * lambda_scale
canvas = ROOT.TCanvas("cCompareRside95", "", 1500, 530)
canvas.Divide(3, 1, 0.006, 0.003)
keep = []
for panel, axis in enumerate(("out", "side", "long"), 1):
    data = fit.Get(f"hDataCk3D_mixedWeighted_proj_{axis}")
    baseline = fit.Get(f"hCATS3D_background_mixedWeighted_proj_{axis}")
    keep += [data, baseline]
    canvas.cd(panel)
    ROOT.gPad.SetLeftMargin(0.15)
    ROOT.gPad.SetRightMargin(0.035)
    ROOT.gPad.SetBottomMargin(0.15)
    ROOT.gPad.SetTopMargin(0.27)
    frame = ROOT.gPad.DrawFrame(-0.04, 0.99, 0.04, 1.22)
    frame.SetTitle(f";q_{{{axis}}} (GeV/c);C(q)")
    frame.GetYaxis().SetTitleOffset(1.15)
    hists = []
    for imodel, (name, array) in enumerate(models.items()):
        hist = ROOT.TH1D(f"hRside_{panel}_{imodel}", "", 4, -0.04, 0.04)
        hist.SetDirectory(0)
        for side in (0, 1):
            sels = [slice(6, 14)] * 3
            sels[panel-1] = slice(10 + 4*side, 14 + 4*side)
            raw = array[tuple(sels)].mean()
            for sign in (-1, 1):
                q = sign * (0.01 + 0.02*side)
                b = hist.FindBin(q)
                hist.SetBinContent(b, raw * baseline.GetBinContent(baseline.FindBin(q)))
            print(axis, name, 0.01+0.02*side, round(hist.GetBinContent(hist.FindBin(0.01+0.02*side)), 6))
        hist.SetLineColor((ROOT.kBlue+1, ROOT.kRed+1, ROOT.kGreen+2)[imodel])
        hist.SetLineStyle(2 if imodel == 0 else 1)
        hist.SetLineWidth(3)
        hist.Draw("HIST SAME")
        hists.append(hist)
    data.SetMarkerStyle(20)
    data.SetMarkerSize(0.72)
    data.SetLineColor(ROOT.kBlack)
    data.Draw("E1 SAME")
    keep.extend(hists)
    if panel == 1:
        legend = ROOT.TLegend(0.51, 0.73 if all_mode else 0.77, 0.96, 0.96)
        legend.SetBorderSize(0)
        legend.SetFillStyle(0)
        legend.SetTextSize(0.030)
        legend.AddEntry(data, "Data", "lep")
        legend.AddEntry(hists[0], "original R", "l")
        legend.AddEntry(hists[1], "all R #times 0.95" if all_mode else "R_{side}=7.60 fm", "l")
        if all_mode:
            legend.AddEntry(hists[2], "all R #times 0.95, peak aligned", "l")
        legend.Draw()
        keep.append(legend)
    label = ROOT.TLatex()
    label.SetNDC()
    label.SetTextSize(0.032)
    label.DrawLatex(0.16, 0.93, "0-10%, 0.2 < m_{T} < 0.4")
    label.SetTextSize(0.025)
    label.DrawLatex(0.16, 0.84, "#omega=0.16; toy q response")
    keep.append(label)

canvas.SaveAs(base + ("toy_Rall95_omega016_qresponse.png" if all_mode
                      else "toy_Rside95_omega016_qresponse.png"))
