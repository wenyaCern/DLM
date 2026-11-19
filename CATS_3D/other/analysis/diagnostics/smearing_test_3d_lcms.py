"""Estimate a 3D C(q) response using the available per-axis true/reco matrices.

This is a factorized cent/mT-integrated diagnostic, not a final 3D response.
"""

from pathlib import Path
import json
import math

import numpy as np
import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
here = Path(__file__).resolve().parents[3]
input_path = here.parent.parent / "3D_pion_CATS/ROOT/lambdaPara/AnalysisResults_lambdaPara_smearing3DLCMS_test.root"
model_path = here / "fit/results/derivedFull_Bminus_cent0to10_anyCPR_fineKernel_derivedLambda_baseline100to200_fit/mt0/qn0/fit_phi0_rinit8p0.root"
output_dir = here / "resolution_test_3d_lcms"
output_dir.mkdir(exist_ok=True)

source = ROOT.TFile.Open(str(input_path))
fit = ROOT.TFile.Open(str(model_path))
if not source or source.IsZombie() or not fit or fit.IsZombie():
    raise RuntimeError("Missing response or model ROOT file")

matrices = []
true_marginals = []
stats = {}
names = ("out", "side", "long")
base = "femto-pair-track-track/TrackTrack/SE/MC/"
for axis in names:
    hist = source.Get(base + f"hTrueQ{axis}VsQ{axis}")
    if not hist or hist.GetNbinsX() != hist.GetNbinsY():
        raise RuntimeError("Unexpected response binning")
    n=hist.GetNbinsX()
    fine_edges=np.array([hist.GetXaxis().GetBinLowEdge(i + 1) for i in range(n + 1)])
    if abs(fine_edges[0]+0.2)>1e-9 or abs(fine_edges[-1]-0.2)>1e-9:
        raise RuntimeError("Unexpected response range")
    counts = np.array([[hist.GetBinContent(i + 1, j + 1) for j in range(n)]
                       for i in range(n)], dtype=float)
    # Split each fine histogram cell by its geometric overlap with the
    # 20-MeV model cells, retaining true/reco counts rather than bin centers.
    coarse_edges=np.linspace(-0.2,0.2,21)
    overlap=np.maximum(0,np.minimum(fine_edges[1:,None],coarse_edges[None,1:])-
                       np.maximum(fine_edges[:-1,None],coarse_edges[None,:-1]))
    overlap/=np.diff(fine_edges)[:,None]
    counts20=overlap.T @ counts @ overlap
    truth = counts20.sum(axis=1)
    matrix = counts20 / np.maximum(truth[:, None], 1.0)
    matrices.append(matrix)
    true_marginals.append(truth / truth.mean())
    central = np.flatnonzero(np.abs((fine_edges[:-1]+fine_edges[1:])/2)<0.04)
    in_range = counts[central].sum()
    same = sum(counts[i, i] for i in central) / in_range
    nearest = sum(counts[i, j] for i in central for j in range(n)
                  if abs(j - i) == 1) / in_range
    stats[axis] = {"entries": hist.GetEntries(),
                   "response_bin_width_MeV": float(1000*np.diff(fine_edges)[0]),
                   "central_same_bin": same,
                   "central_nearest_bin": nearest,
                   "central_farther_than_one_bin": 1 - same - nearest}

hist = fit.Get("hCATS3D_shape")
data_hist = fit.Get("hDataCk3D")
if not hist or any(ax.GetNbins() != 20 for ax in (hist.GetXaxis(), hist.GetYaxis(), hist.GetZaxis())):
    raise RuntimeError("Unexpected model binning")
core = np.array([[[hist.GetBinContent(i + 1, j + 1, k + 1)
                   for k in range(20)] for j in range(20)] for i in range(20)])
data = np.array([[[data_hist.GetBinContent(i + 1, j + 1, k + 1)
                   for k in range(20)] for j in range(20)] for i in range(20)])

def smear(signal, prior, response=matrices):
    # Matrix rows: true bin; columns: reconstructed bin.
    contraction = lambda a: np.einsum("abc,ai,bj,ck->ijk", a, *response, optimize=True)
    return contraction(signal * prior) / np.maximum(contraction(prior), 1e-20)

def gaussian_response(sigma_gev):
    edges=np.linspace(-0.2,0.2,21)
    result=np.empty((20,20))
    offsets=np.linspace(0.025,0.975,40)
    for i in range(20):
        true_positions=edges[i]+offsets*(edges[i+1]-edges[i])
        for j in range(20):
            result[i,j]=sum(0.5*(math.erf((edges[j+1]-t)/(math.sqrt(2)*sigma_gev))-
                                 math.erf((edges[j]-t)/(math.sqrt(2)*sigma_gev)))
                            for t in true_positions)/len(true_positions)
    return result

flat = smear(core, np.ones_like(core))
prior = np.einsum("i,j,k->ijk", *true_marginals)
weighted = smear(core, prior)
toy_cases = {}
for label,sigmas in {"toy_4_4_4_MeV":(4,4,4),
                     "toy_6_4_4_MeV":(6,4,4),
                     "toy_8_6_6_MeV":(8,6,6)}.items():
    toy_cases[label]=smear(core,np.ones_like(core),
                           [gaussian_response(s/1000) for s in sigmas])
cent = (9, 10)
central_cube = np.ix_(cent, cent, cent)
stats["model"] = {
    "central_c_original": float(core[central_cube].mean()),
    "central_c_flat_prior": float(flat[central_cube].mean()),
    "central_c_mc_marginal_prior": float(weighted[central_cube].mean()),
    "central_c_data": float(data[central_cube].mean()),
}
for label,values in toy_cases.items():
    stats["model"]["central_c_"+label]=float(values[central_cube].mean())

canvas = ROOT.TCanvas("cResponse", "True pi-pi factorized MC response", 1800, 480)
canvas.Divide(4, 1)
keep = []
legend_entries = None
for a in range(3):
    def line(arr):
        if a == 0: return arr[:, 9:11, 9:11].mean(axis=(1, 2))
        if a == 1: return arr[9:11, :, 9:11].mean(axis=(0, 2))
        return arr[9:11, 9:11, :].mean(axis=(0, 1))
    pad = canvas.cd(a + 1)
    pad.SetLeftMargin(0.13)
    pad.SetBottomMargin(0.14)
    curves=[]
    for label, values, color in [
        ("Data",line(data),ROOT.kBlack),
        ("Original model",line(core),ROOT.kBlue+1),
        ("Toy core (6,4,4 MeV)",line(toy_cases["toy_6_4_4_MeV"]),ROOT.kGreen+2),
        ("Smeared, flat prior",line(flat),ROOT.kOrange+7),
        ("Smeared, MC prior",line(weighted),ROOT.kViolet+1)]:
        h=ROOT.TH1D(f"axis{a}_{len(curves)}",f";q_{{{names[a]}}} (GeV/c);C(q)",20,-0.2,0.2)
        h.SetDirectory(0)
        for i,value in enumerate(values,1): h.SetBinContent(i,float(value))
        h.SetLineColor(color)
        h.SetMarkerColor(color)
        h.SetLineWidth(2)
        if label=="Data": h.SetMarkerStyle(20); h.SetMarkerSize(0.6)
        if label=="Smeared, MC prior": h.SetLineStyle(2)
        h.GetXaxis().SetRangeUser(-0.1,0.1)
        h.SetMinimum(0.98)
        h.SetMaximum(1.34)
        h.Draw("P" if label=="Data" else "HIST SAME")
        curves.append((h,label))
    if a == 0: legend_entries = curves
    keep.extend(h for h,_ in curves)
canvas.cd(4)
legend=ROOT.TLegend(0.07,0.66,0.95,0.93)
legend.SetBorderSize(0)
legend.SetTextSize(0.045)
for h,label in legend_entries: legend.AddEntry(h,label,"p" if label=="Data" else "l")
legend.Draw()
note=ROOT.TLatex()
note.SetNDC(True)
note.SetTextSize(0.039)
note.DrawLatex(0.07,0.51,"True #pi#pi MC response")
note.DrawLatex(0.07,0.45,"0-80% cent; 0.2 < m_{T} < 0.7")
note.DrawLatex(0.07,0.39,"Factorized 1D q responses")
keep.extend([legend,note])
canvas.SaveAs(str(output_dir / "smearing_shape_estimate.png"))
with open(output_dir / "summary.json", "w") as file:
    json.dump(stats, file, indent=2)
print(json.dumps(stats, indent=2))
