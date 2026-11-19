#!/usr/bin/env python3
import math
import re
from array import array
from pathlib import Path
import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
BASE = Path(__file__).resolve().parents[2]
LOG_BASE = BASE / "fit/nuclearStructureDatafit/results/oo_merged_mT/logs"
OUT_DIR = BASE / "fit/nuclearStructureDatafit/results/oo_merged_mT/summary"
OUT_DIR.mkdir(parents=True, exist_ok=True)
PHI = [0.1, 0.7, 1.3, 1.9, 2.5, 3.1]
COMPONENTS = [(1, "Rout", "R_{out}^{2}"),
              (2, "Rside", "R_{side}^{2}"),
              (3, "Rlong", "R_{long}^{2}")]
QSETS = [(0, "low q_{2} (0-30%)", ROOT.kBlue + 1, 24),
         (2, "high q_{2} (70-100%)", ROOT.kRed + 1, 20)]

def read_component(qn, index, name):
    radii, errors = [], []
    for iphi in range(6):
        text = (LOG_BASE / f"cent1_qn{qn}" / f"phi{iphi}.log").read_text()
        matches = re.findall(
            rf"^\s*{index}\s+{name}\s+([0-9.eE+-]+)\s+([0-9.eE+-]+)",
            text, flags=re.MULTILINE)
        if not matches:
            raise RuntimeError(f"Missing {name} in qn{qn} phi{iphi}")
        value, error = map(float, matches[-1])
        radii.append(value)
        errors.append(error)
    return [r*r for r in radii], [2*r*e for r, e in zip(radii, errors)]

canvas = ROOT.TCanvas("cOOLowHigh", "O-O low/high q2", 1600, 600)
canvas.Divide(3, 1, 0.002, 0.002)
keep = []
summary = []

for ipad, (index, name, ylabel) in enumerate(COMPONENTS, start=1):
    pad = canvas.cd(ipad)
    pad.SetLeftMargin(0.16)
    pad.SetRightMargin(0.035)
    pad.SetBottomMargin(0.17)
    pad.SetTopMargin(0.12)
    pad.SetTicks(1, 1)
    legend = ROOT.TLegend(0.49, 0.72, 0.93, 0.88)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.SetTextSize(0.037)

    graphs = []
    all_y, all_ey = [], []
    for qn, qlabel, color, marker in QSETS:
        y, ey = read_component(qn, index, name)
        all_y.extend(y)
        all_ey.extend(ey)
        graph = ROOT.TGraphErrors(6, array("d", PHI), array("d", y),
                                  array("d", [0.0]*6), array("d", ey))
        graph.SetName(f"g_{name}_q{qn}")
        graph.SetMarkerStyle(marker)
        graph.SetMarkerSize(1.15)
        graph.SetMarkerColor(color)
        graph.SetLineColor(color)
        graphs.append((qn, qlabel, color, graph))

    ymin = min(y-e for y, e in zip(all_y, all_ey))
    ymax = max(y+e for y, e in zip(all_y, all_ey))
    span = max(ymax-ymin, 0.1)
    frame = pad.DrawFrame(-0.2, ymin-0.18*span, 3.4, ymax+0.38*span)
    frame.SetTitle(f";#Delta#varphi = #varphi_{{pair}} - #Psi_{{2}} (rad);{ylabel} (fm^{{2}})")
    frame.GetXaxis().SetTitleSize(0.055)
    frame.GetXaxis().SetLabelSize(0.043)
    frame.GetXaxis().SetTitleOffset(1.25)
    frame.GetYaxis().SetTitleSize(0.057)
    frame.GetYaxis().SetLabelSize(0.043)
    frame.GetYaxis().SetTitleOffset(1.25)
    title = ROOT.TLatex()
    title.SetNDC(True)
    title.SetTextSize(0.050)
    title.DrawLatex(0.18, 0.91, ylabel)
    keep.extend([frame, title, legend])

    for qn, qlabel, color, graph in graphs:
        graph.Draw("P SAME")
        fit = ROOT.TF1(f"fit_{name}_q{qn}", "[0]+[1]*cos(2*x)", -0.2, 3.4)
        fit.SetParameters(sum(graph.GetY()[i] for i in range(6))/6., 0.1)
        fit.SetLineColor(color)
        fit.SetLineWidth(2)
        if qn == 0:
            fit.SetLineStyle(2)
        graph.Fit(fit, "Q0")
        fit.Draw("SAME")
        a, ea = fit.GetParameter(0), fit.GetParError(0)
        b, eb = fit.GetParameter(1), fit.GetParError(1)
        rel = b/a
        erel = math.sqrt((eb/a)**2 + (b*ea/a**2)**2)
        legend.AddEntry(graph, qlabel, "lep")
        summary.append((name, qn, a, ea, b, eb, rel, erel,
                        fit.GetChisquare(), fit.GetNDF()))
        keep.extend([graph, fit])
    legend.Draw()

png = OUT_DIR / "OO_cent10to30_low_vs_high_q2_R2_vs_phi.png"
pdf = OUT_DIR / "OO_cent10to30_low_vs_high_q2_R2_vs_phi.pdf"
canvas.SaveAs(str(png))
canvas.SaveAs(str(pdf))

csv = OUT_DIR / "OO_cent10to30_low_vs_high_q2_R2_fourier.csv"
with csv.open("w") as out:
    out.write("component,qn,A,Aerr,B,Berr,B_over_A,B_over_A_err,chi2,ndf\n")
    for row in summary:
        out.write(",".join(map(str, row)) + "\n")
print(png)
