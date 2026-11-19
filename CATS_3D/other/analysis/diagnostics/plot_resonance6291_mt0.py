import csv

import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
base = "derivedFull_Bminus_cent0to10_anyCPR_resonance6291_q10_r20_weak1557_fit/mt0/qn0/"
reso = ROOT.TFile.Open(base + "fit_phi0.root")
gauss = ROOT.TFile.Open(base + "gaussian_control_phi0.root")
with open(base + "r20_vs_r40_selected_q.csv", newline="") as stream:
    selected = {(row["axis"], float(row["qcenter_GeV"])): float(row["modelC"])
                for row in csv.DictReader(stream) if float(row["rmax_fm"]) == 40.0}
with open(base + "gaussian_r20_vs_r40_selected_q.csv", newline="") as stream:
    selected_gauss = {(row["axis"], float(row["qcenter_GeV"])): float(row["modelC"])
                      for row in csv.DictReader(stream) if float(row["rmax_fm"]) == 40.0}

canvas = ROOT.TCanvas("cResonance6291", "", 1620, 650)
canvas.Divide(3, 1, 0.003, 0.003)
keep = []
for panel, axis in enumerate(("out", "side", "long"), 1):
    canvas.cd(panel)
    ROOT.gPad.SetLeftMargin(0.15)
    ROOT.gPad.SetRightMargin(0.03)
    ROOT.gPad.SetBottomMargin(0.15)
    ROOT.gPad.SetTopMargin(0.28)
    frame = ROOT.gPad.DrawFrame(-0.10, 0.92, 0.10, 1.35)
    frame.SetTitle(f";q_{{{axis}}} (GeV/c);C(q)")
    frame.GetYaxis().SetTitleOffset(1.25)
    data = reso.Get(f"hDataCk3D_mixedWeighted_proj_{axis}")
    model_r20 = reso.Get(f"hCATS3D_best_mixedWeighted_proj_{axis}")
    model_gauss = gauss.Get(f"hCATS3D_best_mixedWeighted_proj_{axis}")
    baseline = reso.Get(f"hCATS3D_background_mixedWeighted_proj_{axis}")
    model_gauss.SetLineColor(ROOT.kBlue+1)
    model_gauss.SetLineWidth(3)
    model_gauss.Draw("HIST SAME")
    model_r20.SetLineColor(ROOT.kRed+1)
    model_r20.SetLineWidth(3)
    model_r20.Draw("HIST SAME")
    points = ROOT.TGraph(4)
    for i, q in enumerate((-0.03, -0.01, 0.01, 0.03)):
        b = baseline.GetBinContent(baseline.FindBin(q))
        points.SetPoint(i, q, selected[(axis, abs(q))] * b)
    points.SetMarkerStyle(29)
    points.SetMarkerSize(1.6)
    points.SetMarkerColor(ROOT.kGreen+2)
    points.Draw("P SAME")
    data.SetMarkerStyle(20)
    data.SetMarkerSize(0.8)
    data.SetLineColor(ROOT.kBlack)
    data.Draw("E1 SAME")
    label = ROOT.TLatex()
    label.SetNDC()
    label.SetTextSize(0.033)
    label.DrawLatex(0.16, 0.95, "0-10%; 0.2 < m_{T} < 0.4; qn0, #phi0")
    keep.extend((data, model_r20, model_gauss, baseline, points, label))
    if panel == 1:
        legend = ROOT.TLegend(0.12, 0.75, 0.96, 0.92)
        legend.SetBorderSize(0)
        legend.SetFillStyle(0)
        legend.SetTextSize(0.030)
        legend.AddEntry(data, "Data", "lep")
        legend.AddEntry(model_gauss, "Gaussian, r_{max}=20 fm", "l")
        legend.AddEntry(model_r20, "Resonance 62.91%, r_{max}=20 fm", "l")
        legend.AddEntry(points, "Resonance, r_{max}=40 fm (selected q)", "p")
        legend.Draw()
        keep.append(legend)

canvas.SaveAs(base + "resonance6291_vs_gaussian_slices.png")

compare = ROOT.TCanvas("cR40Selected", "", 1620, 630)
compare.Divide(3, 1, 0.003, 0.003)
keep_r40 = []
for panel, axis in enumerate(("out", "side", "long"), 1):
    compare.cd(panel)
    ROOT.gPad.SetLeftMargin(0.15)
    ROOT.gPad.SetRightMargin(0.03)
    ROOT.gPad.SetBottomMargin(0.15)
    ROOT.gPad.SetTopMargin(0.26)
    frame = ROOT.gPad.DrawFrame(-0.045, 1.02, 0.045, 1.255)
    frame.SetTitle(f";q_{{{axis}}} (GeV/c);C(q)")
    frame.GetYaxis().SetTitleOffset(1.25)
    data = reso.Get(f"hDataCk3D_mixedWeighted_proj_{axis}")
    baseline = reso.Get(f"hCATS3D_background_mixedWeighted_proj_{axis}")
    graphs = []
    for values, color, style in ((selected_gauss, ROOT.kBlue+1, 25),
                                 (selected, ROOT.kGreen+2, 29)):
        graph = ROOT.TGraph(4)
        for i, q in enumerate((-0.03, -0.01, 0.01, 0.03)):
            b = baseline.GetBinContent(baseline.FindBin(q))
            graph.SetPoint(i, q, values[(axis, abs(q))] * b)
        graph.SetMarkerColor(color)
        graph.SetMarkerStyle(style)
        graph.SetMarkerSize(1.7)
        graph.Draw("P SAME")
        graphs.append(graph)
    data.SetMarkerStyle(20)
    data.SetMarkerSize(0.95)
    data.SetLineColor(ROOT.kBlack)
    data.Draw("E1 SAME")
    label = ROOT.TLatex()
    label.SetNDC()
    label.SetTextSize(0.032)
    label.DrawLatex(0.16, 0.94, "0-10%; 0.2 < m_{T} < 0.4; qn0, #phi0")
    keep_r40.extend((data, baseline, *graphs, label))
    if panel == 1:
        legend = ROOT.TLegend(0.12, 0.76, 0.95, 0.91)
        legend.SetBorderSize(0)
        legend.SetFillStyle(0)
        legend.SetTextSize(0.034)
        legend.AddEntry(data, "Data", "lep")
        legend.AddEntry(graphs[0], "Gaussian, r_{max}=40 fm", "p")
        legend.AddEntry(graphs[1], "Resonance 62.91%, r_{max}=40 fm", "p")
        legend.Draw()
        keep_r40.append(legend)

compare.SaveAs(base + "gaussian_vs_resonance_r40_selected_q.png")
