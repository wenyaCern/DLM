"""Fixed-parameter Gaussian radius/lambda scan with the existing toy q response."""
import csv
import math

import numpy as np
import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
base = "fit/results/derivedFull_Bminus_cent0to10_anyCPR_fineKernel_derivedLambda_baseline100to200_fit/mt0/qn0/"
fit = ROOT.TFile.Open(base + "fit_phi0_lambda083_rinit8p0.root")
old_lambda = 0.83 * 0.86**2
weak_fraction = 0.1557
q_mapping = np.array([9-i if i < 10 else i-10 for i in range(20)])


def blur_axis(array, sigma_mev, axis):
    sigma = sigma_mev / 5.0
    radius = math.ceil(4 * sigma)
    offsets = np.arange(-radius, radius + 1)
    weights = np.exp(-0.5 * (offsets / sigma)**2)
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


def load_and_smear(scale):
    tag = "" if scale == 1.0 else f"_Rall{round(scale * 100):03d}"
    octant = np.empty((10, 10, 10), dtype=float)
    with open(base + f"cats3d_r40_lambda083_q5{tag}_octant.csv", newline="") as stream:
        for row in csv.DictReader(stream):
            octant[int(row["i_out"]), int(row["i_side"]), int(row["i_long"])] = float(row["C"])
    model = octant[np.ix_(q_mapping, q_mapping, q_mapping)]
    for axis, sigma in enumerate((8.0, 2.6, 2.6)):
        model = blur_axis(model, sigma, axis)
    return model


def bin_mean(model, axis, shoulder):
    selections = [slice(6, 14)] * 3
    selections[axis] = slice(10 + 4*shoulder, 14 + 4*shoulder)
    return model[tuple(selections)].mean()


rows = []
for scale in (0.90, 0.95, 1.0, 1.05, 1.10):
    model = load_and_smear(scale)
    data_out = fit.Get("hDataCk3D_mixedWeighted_proj_out")
    base_out = fit.Get("hCATS3D_background_mixedWeighted_proj_out")
    target = data_out.GetBinContent(data_out.FindBin(0.01))
    bpeak = base_out.GetBinContent(base_out.FindBin(0.01))
    matched_lambda = (target / bpeak - 1) * old_lambda / ((1-weak_fraction)**2 * (bin_mean(model, 0, 0)-1))
    for lambda_gen in (0.83, 0.86, 0.90, matched_lambda):
        for axis, label in enumerate(("out", "side", "long")):
            data = fit.Get(f"hDataCk3D_mixedWeighted_proj_{label}")
            baseline = fit.Get(f"hCATS3D_background_mixedWeighted_proj_{label}")
            for shoulder in (0, 1):
                q = 0.01 + 0.02*shoulder
                b = baseline.GetBinContent(baseline.FindBin(q))
                prediction = b * (1 + lambda_gen*(1-weak_fraction)**2 / old_lambda * (bin_mean(model, axis, shoulder)-1))
                rows.append((scale, "peak_matched" if lambda_gen == matched_lambda else "fixed", lambda_gen,
                             label, q, prediction, data.GetBinContent(data.FindBin(q))))

output = base + "gaussian_R_scan_omega1557_toy_response.csv"
with open(output, "w", newline="") as stream:
    writer = csv.writer(stream)
    writer.writerow(("R_scale", "lambda_mode", "lambda_gen", "axis", "q_GeV", "model_C", "data_C"))
    writer.writerows(rows)
for row in rows:
    if row[3] == "out":
        print("R=%.2f mode=%s lambda=%.5f q=%.2f model=%.6f data=%.6f" % (row[0], row[1], row[2], row[4], row[5], row[6]))

canvas = ROOT.TCanvas("scan", "Gaussian radius scan", 1500, 500)
canvas.Divide(3, 1)
keep = []
for pad_index, label in enumerate(("out", "side", "long"), 1):
    canvas.cd(pad_index)
    ROOT.gPad.SetLeftMargin(0.13)
    ROOT.gPad.SetBottomMargin(0.13)
    frame = ROOT.TH1F(f"frame_{label}", f"q_{{{label}}};q_{{{label}}} (GeV/c);C(q)", 4, -0.04, 0.04)
    frame.SetMinimum(1.0)
    frame.SetMaximum(1.22)
    frame.Draw("AXIS")
    keep.append(frame)
    legend = ROOT.TLegend(0.43, 0.60, 0.96, 0.94)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    for scale, color in ((0.90, ROOT.kMagenta+1), (0.95, ROOT.kGreen+2), (1.0, ROOT.kBlue+1), (1.05, ROOT.kOrange+7), (1.10, ROOT.kRed+1)):
        matching = [r for r in rows if r[0] == scale and r[1] == "peak_matched" and r[3] == label]
        peak = next(r[5] for r in matching if r[4] == 0.01)
        shoulder = next(r[5] for r in matching if r[4] == 0.03)
        graph = ROOT.TGraph(4)
        for i, (q, value) in enumerate(zip((-0.03, -0.01, 0.01, 0.03), (shoulder, peak, peak, shoulder))):
            graph.SetPoint(i, q, value)
        graph.SetLineColor(color)
        graph.SetMarkerColor(color)
        graph.SetMarkerStyle(20)
        graph.SetLineWidth(2)
        graph.Draw("LP SAME")
        legend.AddEntry(graph, f"R x {scale:.2f}, #lambda={matching[0][2]:.3f}", "lp")
        keep.append(graph)
    data_hist = fit.Get(f"hDataCk3D_mixedWeighted_proj_{label}")
    data_graph = ROOT.TGraphErrors(4)
    for i, q in enumerate((-0.03, -0.01, 0.01, 0.03)):
        b = data_hist.FindBin(q)
        data_graph.SetPoint(i, q, data_hist.GetBinContent(b))
        data_graph.SetPointError(i, 0, data_hist.GetBinError(b))
    data_graph.SetMarkerStyle(20)
    data_graph.SetMarkerSize(1.2)
    data_graph.Draw("P SAME")
    legend.AddEntry(data_graph, "Data", "p")
    legend.Draw()
    keep.extend((data_graph, legend))
canvas.SaveAs(base + "gaussian_R_peak_matched_toy_response.png")
