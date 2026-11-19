"""Compare the qlong~0 out-side slice before/after the diagonal baseline term."""

from pathlib import Path
import ROOT


ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetPalette(ROOT.kBird)

directory = Path(__file__).resolve().parents[3] / (
    "fit/results/derivedFull_Bminus_cent0to10_anyCPR_adaptiveKernel_smear3d_weak155_fit"
    "/mt0/qn0"
)
files = [
    ROOT.TFile.Open(str(directory / "fit_phi0_PIDresponse.root")),
    ROOT.TFile.Open(str(directory / "fit_phi0_PIDresponse_diagonalBand_fixed_kernelFitBoxOnly.root")),
]
data = files[0].Get("hDataCk3D")
old_model = files[0].Get("hCATS3D_best")
band_model = files[1].Get("hCATS3D_best")


def project(source, name):
    result = ROOT.TH2D(name, ";q_{out} (GeV/c);q_{side} (GeV/c)",
                       20, -0.2, 0.2, 20, -0.2, 0.2)
    for ix in range(1, 21):
        for iy in range(1, 21):
            num = den = 0.0
            for iz in (10, 11):  # |qlong| < 0.02 GeV/c
                value = source.GetBinContent(ix, iy, iz)
                error = data.GetBinError(ix, iy, iz)
                if error <= 0:
                    continue
                weight = 1.0 / (error * error)
                num += weight * value
                den += weight
            if den:
                result.SetBinContent(ix, iy, num / den)
    return result


h_data = project(data, "hData2D")
h_old = project(old_model, "hOld2D")
h_new = project(band_model, "hBand2D")
h_old_res = h_data.Clone("hOldResidual2D")
h_old_res.Divide(h_old)
h_new_res = h_data.Clone("hBandResidual2D")
h_new_res.Divide(h_new)
for ix in range(1, 21):
    for iy in range(1, 21):
        for h in (h_old_res, h_new_res):
            if h.GetBinContent(ix, iy):
                h.SetBinContent(ix, iy, h.GetBinContent(ix, iy) - 1)

canvas = ROOT.TCanvas("band_canvas", "band_canvas", 1500, 1080)
canvas.Divide(2, 2)
labels = [
    (h_data, "Data", 0.97, 1.20),
    (h_new, "Fit region: kernel; outside: baseline", 0.97, 1.20),
    (h_old_res, "Data / fit - 1; no band", -0.025, 0.025),
    (h_new_res, "Data / fit - 1; band fixed from sideband", -0.025, 0.025),
]
keep = []
for pad, (hist, title, low, high) in enumerate(labels, 1):
    canvas.cd(pad)
    ROOT.gPad.SetLeftMargin(0.12)
    ROOT.gPad.SetRightMargin(0.16)
    ROOT.gPad.SetTopMargin(0.11)
    hist.SetTitle("")
    hist.SetMinimum(low)
    hist.SetMaximum(high)
    hist.GetXaxis().SetTitleSize(0.045)
    hist.GetYaxis().SetTitleSize(0.045)
    hist.Draw("COLZ")
    line = ROOT.TLine(-0.2, -0.03, 0.2, 0.03)
    line.SetLineStyle(2)
    line.SetLineWidth(2)
    line.Draw()
    text = ROOT.TLatex()
    text.SetNDC(True)
    text.SetTextSize(0.045)
    text.DrawLatex(0.12, 0.92, title)
    keep.extend([line, text])

out = directory / "diagonal_band_baseline_2D_comparison_fitBoxOnly.png"
canvas.SaveAs(str(out))
print(out)
