"""Visualize fitted q-bin averages with PCHIP curves and a radius-covariance band.

The smooth curves interpolate saved bin averages; they are not extra theory
evaluations between model grid points. The band propagates the three fitted
radius errors and correlations, holding lambda and baseline fixed.
"""
import argparse
import math
from pathlib import Path

import numpy as np
import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)


def pchip(x, y, xx):
    """Shape-preserving cubic Hermite interpolation on strictly increasing x."""
    h = np.diff(x)
    delta = np.diff(y) / h
    d = np.zeros_like(y)
    for k in range(1, len(y) - 1):
        if delta[k - 1] * delta[k] > 0:
            w1 = 2 * h[k] + h[k - 1]
            w2 = h[k] + 2 * h[k - 1]
            d[k] = (w1 + w2) / (w1 / delta[k - 1] + w2 / delta[k])
    for k, side in ((0, 0), (len(y) - 1, -1)):
        if k == 0:
            slope = ((2 * h[0] + h[1]) * delta[0] - h[0] * delta[1]) / (h[0] + h[1])
            first, second = delta[0], delta[1]
        else:
            slope = ((2 * h[-1] + h[-2]) * delta[-1] - h[-1] * delta[-2]) / (h[-1] + h[-2])
            first, second = delta[-1], delta[-2]
        if np.sign(slope) != np.sign(first):
            slope = 0
        elif np.sign(first) != np.sign(second) and abs(slope) > 3 * abs(first):
            slope = 3 * first
        d[k] = slope
    idx = np.clip(np.searchsorted(x, xx, side="right") - 1, 0, len(h) - 1)
    t = (xx - x[idx]) / h[idx]
    return ((2*t**3 - 3*t**2 + 1) * y[idx]
            + (t**3 - 2*t**2 + t) * h[idx] * d[idx]
            + (-2*t**3 + 3*t**2) * y[idx + 1]
            + (t**3 - t**2) * h[idx] * d[idx + 1])


def make_graph(x, y, name):
    graph = ROOT.TGraph(len(x))
    graph.SetName(name)
    for i, (a, b) in enumerate(zip(x, y)):
        graph.SetPoint(i, float(a), float(b))
    return graph


parser = argparse.ArgumentParser()
parser.add_argument("fit_root")
parser.add_argument("output_png")
args = parser.parse_args()

file = ROOT.TFile.Open(args.fit_root)
if not file or file.IsZombie():
    raise RuntimeError(f"Cannot open {args.fit_root}")
cov_hist = file.Get("hRadiusCovariance")
if not cov_hist:
    raise RuntimeError("Fit ROOT file does not contain radius covariance")
cov = np.array([[cov_hist.GetBinContent(a + 1, b + 1) for b in range(3)]
                for a in range(3)])
if np.linalg.eigvalsh(cov).min() < -1e-10:
    raise RuntimeError("Radius covariance is not positive semidefinite")

canvas = ROOT.TCanvas("smoothRadiusBand", "", 1650, 560)
canvas.Divide(3, 1, 0.006, 0.003)
objects = []
max_sigmas = {}
uncertainty_curves = {}
for panel, axis in enumerate(("out", "side", "long"), 1):
    data = file.Get(f"hDataCk3D_mixedWeighted_proj_{axis}")
    fit = file.Get(f"hCATS3D_best_mixedWeighted_proj_{axis}")
    baseline = file.Get(f"hCATS3D_background_mixedWeighted_proj_{axis}")
    derivatives = [file.Get(f"hCATS3D_d{name}_mixedWeighted_proj_{axis}")
                   for name in ("Rout", "Rside", "Rlong")]
    if not all([data, fit, baseline, *derivatives]):
        raise RuntimeError(f"Missing {axis} projection or derivative")
    x = np.array([fit.GetBinCenter(i) for i in range(1, fit.GetNbinsX() + 1)])
    y = np.array([fit.GetBinContent(i) for i in range(1, fit.GetNbinsX() + 1)])
    grad = np.array([[h.GetBinContent(i) for h in derivatives]
                     for i in range(1, fit.GetNbinsX() + 1)])
    sigma = np.sqrt(np.maximum(0, np.einsum("ia,ab,ib->i", grad, cov, grad)))
    max_sigmas[axis] = float(sigma.max())
    xx = np.linspace(x[0], x[-1], 801)
    uncertainty_curves[axis] = (xx, pchip(x, sigma, xx))
    center = pchip(x, y, xx)
    upper = pchip(x, y + sigma, xx)
    lower = pchip(x, y - sigma, xx)
    polygon = make_graph(np.r_[xx, xx[::-1]],
                         np.r_[upper, lower[::-1]], f"band_{axis}")
    polygon.SetFillColorAlpha(ROOT.kCyan + 1, 0.30)
    polygon.SetLineColorAlpha(ROOT.kCyan + 1, 0.0)
    top = make_graph(xx, upper, f"upper_{axis}")
    bottom = make_graph(xx, lower, f"lower_{axis}")
    middle = make_graph(xx, center, f"middle_{axis}")
    for boundary in (top, bottom):
        boundary.SetLineColor(ROOT.kBlue + 1)
        boundary.SetLineWidth(2)
    middle.SetLineColor(ROOT.kMagenta + 1)
    middle.SetLineWidth(3)
    canvas.cd(panel)
    ROOT.gPad.SetLeftMargin(0.15)
    ROOT.gPad.SetRightMargin(0.035)
    ROOT.gPad.SetBottomMargin(0.15)
    ROOT.gPad.SetTopMargin(0.15)
    frame = ROOT.gPad.DrawFrame(-0.2, 0.90, 0.2, 1.30)
    frame.SetTitle(f";q_{{{axis}}} (GeV/c);C(q)")
    frame.GetXaxis().SetTitleSize(0.052)
    frame.GetYaxis().SetTitleSize(0.052)
    frame.GetXaxis().SetLabelSize(0.042)
    frame.GetYaxis().SetLabelSize(0.042)
    frame.GetYaxis().SetTitleOffset(1.22)
    baseline.SetLineColor(ROOT.kGray + 2)
    baseline.SetLineStyle(7)
    baseline.SetLineWidth(2)
    data.SetMarkerStyle(20)
    data.SetMarkerSize(0.65)
    data.SetLineColor(ROOT.kBlack)
    baseline.Draw("HIST SAME")
    polygon.Draw("F SAME")
    top.Draw("L SAME")
    bottom.Draw("L SAME")
    middle.Draw("L SAME")
    data.Draw("E1 SAME")
    label = ROOT.TLatex()
    label.SetNDC(True)
    label.SetTextSize(0.032)
    label.DrawLatex(0.16, 0.91, "0-10%, 0.2 < m_{T} < 0.4; q_{2} bin 0, #phi bin 0")
    label.SetTextSize(0.027)
    label.DrawLatex(0.16, 0.86, "other two |q| < 0.02 GeV/c")
    label.DrawLatex(0.16, 0.81, f"max radius 1#sigma = {1000*sigma.max():.2f} #times 10^{{-3}}")
    if panel == 1:
        legend = ROOT.TLegend(0.61, 0.59, 0.95, 0.79)
        legend.SetBorderSize(0)
        legend.SetFillStyle(0)
        legend.SetTextSize(0.027)
        legend.AddEntry(data, "Data", "lep")
        legend.AddEntry(middle, "interpolated fit", "l")
        legend.AddEntry(polygon, "radius fit #pm1#sigma", "f")
        legend.AddEntry(baseline, "outer-box baseline", "l")
        legend.Draw()
        objects.append(legend)
    objects.extend([polygon, top, bottom, middle, label, data, baseline, frame])

Path(args.output_png).parent.mkdir(parents=True, exist_ok=True)
canvas.SaveAs(args.output_png)
detail = ROOT.TCanvas("radiusBandDetail", "", 1650, 420)
detail.Divide(3, 1, 0.006, 0.003)
detail_objects = []
limit = 1.25 * max(max_sigmas.values())
for panel, axis in enumerate(("out", "side", "long"), 1):
    detail.cd(panel)
    ROOT.gPad.SetLeftMargin(0.15)
    ROOT.gPad.SetRightMargin(0.035)
    ROOT.gPad.SetBottomMargin(0.16)
    ROOT.gPad.SetTopMargin(0.10)
    frame = ROOT.gPad.DrawFrame(-0.2, -limit, 0.2, limit)
    frame.SetTitle(f";q_{{{axis}}} (GeV/c);C_{{fit}} - C_{{center}}")
    frame.GetYaxis().SetTitleOffset(1.30)
    xx, smooth_sigma = uncertainty_curves[axis]
    top = make_graph(xx, smooth_sigma, f"detailUpper_{axis}")
    bottom = make_graph(xx, -smooth_sigma, f"detailLower_{axis}")
    zero = make_graph(np.array([-0.2, 0.2]), np.array([0.0, 0.0]), f"detailCenter_{axis}")
    for line in (top, bottom):
        line.SetLineColor(ROOT.kBlue + 1)
        line.SetLineWidth(2)
        line.Draw("L SAME")
    zero.SetLineColor(ROOT.kMagenta + 1)
    zero.SetLineWidth(2)
    zero.Draw("L SAME")
    detail_objects.extend([frame, top, bottom, zero])
detail_path = str(Path(args.output_png).with_name(Path(args.output_png).stem + "_uncertainty_zoom.png"))
detail.SaveAs(detail_path)
print("max 1-sigma by axis:", max_sigmas)
