#!/usr/bin/env python3
import math
import re
import sys
from array import array
from pathlib import Path
import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
BASE = Path(__file__).resolve().parents[2]
QN = int(sys.argv[1]) if len(sys.argv) > 1 else 2
if QN not in (0, 1, 2):
    raise ValueError("QN must be 0, 1, or 2")
QN_LABEL = ["low", "mid", "high"][QN]
LOG_DIR = BASE / f"fit/nuclearStructureDatafit/results/oo_merged_mT/logs/cent1_qn{QN}"
OUT_DIR = BASE / "fit/nuclearStructureDatafit/results/oo_merged_mT/summary"
OUT_DIR.mkdir(parents=True, exist_ok=True)
phi = [0.1, 0.7, 1.3, 1.9, 2.5, 3.1]
names = ["Rout", "Rside", "Rlong"]
axis_labels = ["R_{out}^{2}", "R_{side}^{2}", "R_{long}^{2}"]
values = {name: [] for name in names}
errors = {name: [] for name in names}

for iphi in range(6):
    text = (LOG_DIR / f"phi{iphi}.log").read_text()
    for index, name in enumerate(names, start=1):
        matches = re.findall(
            rf"^\s*{index}\s+{name}\s+([0-9.eE+-]+)\s+([0-9.eE+-]+)",
            text, flags=re.MULTILINE)
        if not matches:
            raise RuntimeError(f"Missing final {name} error in phi{iphi}.log")
        value, error = map(float, matches[-1])
        values[name].append(value)
        errors[name].append(error)

canvas = ROOT.TCanvas("cOOQ2", f"O-O {QN_LABEL} q2", 1450, 470)
canvas.Divide(3, 1)
objects, rows = [], []
for ipad, (name, axis_label) in enumerate(zip(names, axis_labels), start=1):
    canvas.cd(ipad)
    ROOT.gPad.SetTicks(1, 1)
    y = [r*r for r in values[name]]
    ey = [2*r*e for r, e in zip(values[name], errors[name])]
    graph = ROOT.TGraphErrors(6, array("d", phi), array("d", y),
                              array("d", [0.0]*6), array("d", ey))
    graph.SetName(f"g{name}2")
    graph.SetTitle(f"{name};#Delta#varphi=#varphi_{{pair}}-#Psi_{{2}} (rad);"
                   f"{axis_label} (fm^{{2}})")
    graph.SetMarkerStyle(20)
    graph.Draw("AP")
    fit = ROOT.TF1(f"f{name}2", "[0]+[1]*cos(2*x)", -0.2, 3.4)
    fit.SetParameters(sum(y)/len(y), 0.1)
    fit.SetLineColor(ROOT.kBlue + 1)
    graph.Fit(fit, "Q")
    a, ea = fit.GetParameter(0), fit.GetParError(0)
    b, eb = fit.GetParameter(1), fit.GetParError(1)
    relative = b/a
    relative_error = math.sqrt((eb/a)**2 + (b*ea/a**2)**2)
    note = ROOT.TLatex()
    note.SetNDC(True)
    note.SetTextSize(0.038)
    note.DrawLatex(0.14, 0.86, f"B/A = {relative:+.4f} #pm {relative_error:.4f}")
    note.DrawLatex(0.14, 0.80,
                   f"#chi^{{2}}/ndf = {fit.GetChisquare():.1f}/{fit.GetNDF()}")
    rows.append((name, a, ea, b, eb, relative, relative_error,
                 fit.GetChisquare(), fit.GetNDF()))
    objects.extend([graph, fit, note])

png = OUT_DIR / f"OO_cent10to30_{QN_LABEL}q2_R2_vs_phi.png"
pdf = OUT_DIR / f"OO_cent10to30_{QN_LABEL}q2_R2_vs_phi.pdf"
canvas.SaveAs(str(png))
canvas.SaveAs(str(pdf))
csv = OUT_DIR / f"OO_cent10to30_{QN_LABEL}q2_R2_fourier.csv"
with csv.open("w") as out:
    out.write("component,A,Aerr,B,Berr,B_over_A,B_over_A_err,chi2,ndf\n")
    for row in rows:
        out.write(",".join(map(str, row)) + "\n")
for row in rows:
    print(f"{row[0]}: A={row[1]:.5f}+/-{row[2]:.5f}, "
          f"B={row[3]:+.5f}+/-{row[4]:.5f}, "
          f"B/A={row[5]:+.4f}+/-{row[6]:.4f}, "
          f"significance={row[3]/row[4]:+.2f} sigma, "
          f"chi2/ndf={row[7]:.2f}/{row[8]}")
