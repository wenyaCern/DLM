"""Exploratory Gaussian q smearing of a fixed-R CATS3D model; not an MC response correction."""
import csv
import math
import sys

import numpy as np
import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
omega = float(sys.argv[1]) if len(sys.argv) > 1 else 0.14
if not 0 <= omega < 1:
    raise ValueError("omega fraction must be in [0,1)")
suffix = "" if abs(omega - 0.14) < 1e-9 else f"_omega{omega:.2f}".replace(".", "")

base = "fit/results/derivedFull_Bminus_cent0to10_anyCPR_fineKernel_derivedLambda_baseline100to200_fit/mt0/qn0/"
fit = ROOT.TFile.Open(base + "fit_phi0_lambda083_rinit8p0.root")
if not fit or fit.IsZombie():
    raise RuntimeError("Cannot open fit ROOT file")

octant = np.empty((10, 10, 10), dtype=float)
with open(base + "cats3d_r40_lambda083_q5_octant.csv", newline="") as stream:
    for row in csv.DictReader(stream):
        octant[int(row["i_out"]), int(row["i_side"]), int(row["i_long"])] = float(row["C"])
mapping = np.array([9-i if i < 10 else i-10 for i in range(20)])
model = octant[np.ix_(mapping, mapping, mapping)]
model = 1.0 + (model - 1.0) * ((1.0 - omega) / 0.86) ** 2


def blur_axis(array, sigma_mev, axis):
    if sigma_mev == 0:
        return array.copy()
    sigma = sigma_mev / 5.0
    radius = max(1, math.ceil(4 * sigma))
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


def smear(array, sigma_out):
    for axis, sigma in enumerate((sigma_out, 2.6, 2.6)):
        array = blur_axis(array, sigma, axis)
    return array


def coarse_slice(array, axis, side):
    # side=0 or 1 corresponds to q_axis = 0–20 or 20–40 MeV/c.
    selections = [slice(6, 14)] * 3  # other two |q| < 20 MeV/c
    selections[axis] = slice(10 + 4 * side, 14 + 4 * side)
    return array[tuple(selections)].mean()


results = {0: model}
for sigma in (4, 6, 8):
    results[sigma] = smear(model, sigma)

canvas = ROOT.TCanvas("cToyQResolution", "", 1500, 530)
canvas.Divide(3, 1, 0.006, 0.003)
keep = []
colors = {0: ROOT.kBlue + 1, 4: ROOT.kGreen + 2,
          6: ROOT.kOrange + 7, 8: ROOT.kRed + 1}
lines = []
for panel, axis in enumerate(("out", "side", "long"), 1):
    data = fit.Get(f"hDataCk3D_mixedWeighted_proj_{axis}")
    baseline = fit.Get(f"hCATS3D_background_mixedWeighted_proj_{axis}")
    keep += [data, baseline]
    canvas.cd(panel)
    ROOT.gPad.SetLeftMargin(0.15)
    ROOT.gPad.SetRightMargin(0.035)
    ROOT.gPad.SetBottomMargin(0.15)
    ROOT.gPad.SetTopMargin(0.28)
    frame = ROOT.gPad.DrawFrame(-0.04, 0.99, 0.04, 1.25)
    frame.SetTitle(f";q_{{{axis}}} (GeV/c);C(q)")
    frame.GetYaxis().SetTitleOffset(1.15)
    for sigma, array in results.items():
        hist = ROOT.TH1D(f"hToy_{axis}_{sigma}", "", 4, -0.04, 0.04)
        hist.SetDirectory(0)
        for side in (0, 1):
            raw = coarse_slice(array, panel - 1, side)
            for sign in (-1, 1):
                q = sign * (0.01 + 0.02 * side)
                b = hist.FindBin(q)
                hist.SetBinContent(b, raw * baseline.GetBinContent(baseline.FindBin(q)))
                if sign > 0:
                    data_y = data.GetBinContent(data.FindBin(q))
                    lines.append((axis, sigma, q, hist.GetBinContent(b), data_y))
        hist.SetLineColor(colors[sigma])
        hist.SetLineWidth(3 if sigma == 0 else 2)
        if sigma == 0:
            hist.SetLineStyle(2)
        hist.Draw("HIST SAME")
        keep.append(hist)
    data.SetMarkerStyle(20)
    data.SetMarkerSize(0.72)
    data.SetLineColor(ROOT.kBlack)
    data.Draw("E1 SAME")
    if panel == 1:
        legend = ROOT.TLegend(0.54, 0.75, 0.96, 0.96)
        legend.SetBorderSize(0)
        legend.SetFillStyle(0)
        legend.SetTextSize(0.029)
        legend.AddEntry(data, "Data", "lep")
        for sigma in (0, 4, 6, 8):
            legend.AddEntry(next(h for h in keep if h.GetName() == f"hToy_{axis}_{sigma}"),
                            "no smearing" if sigma == 0 else f"#sigma_{{qout}}={sigma} MeV/c", "l")
        legend.Draw()
        keep.append(legend)
    label = ROOT.TLatex()
    label.SetNDC()
    label.SetTextSize(0.033)
    label.DrawLatex(0.16, 0.94, f"0-10%, 0.2 < m_{{T}} < 0.4; #omega={omega:.2f}")
    label.SetTextSize(0.025)
    label.DrawLatex(0.16, 0.85, "toy Gaussian response, fixed R")
    keep.append(label)

canvas.SaveAs(base + f"toy_q_resolution_lambda083_r40{suffix}.png")
with open(base + f"toy_q_resolution_lambda083_r40{suffix}.csv", "w", newline="") as stream:
    writer = csv.writer(stream)
    writer.writerow(("axis", "sigma_qout_MeV", "q_center_GeV", "model", "data"))
    writer.writerows(lines)
for line in lines:
    print(*line)
