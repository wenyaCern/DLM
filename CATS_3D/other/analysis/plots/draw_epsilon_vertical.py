"""Stacked mT panels, observed and resolution-corrected eccentricity.
Uses saved HESSE Fourier covariance propagation. Resolution errors are the
previous independent-correlation approximation; epsilon/resolution covariance
is unavailable and set to zero. Same event-weighted resolution for both mT bins.
"""
from pathlib import Path
import json
import csv
import math
import argparse
import ROOT
parser=argparse.ArgumentParser()
parser.add_argument('--xaxis',choices=['percentile','v2'],default='percentile')
parser.add_argument('--bootstrap-errors',action='store_true')
args=parser.parse_args()
# User-provided cent10-30 values, ordered qn0/1/2. Same event selection
# coordinates are used for both mT intervals and both plane measurements.
v2=[0.0703302,0.0720694,0.0760253]
v2err=[6.43559e-05,5.50767e-05,6.09417e-05]
x=v2 if args.xaxis=='v2' else [15,50,85]
ex=v2err if args.xaxis=='v2' else [15,20,15]
xlo,xhi=(0.0695,0.077) if args.xaxis=='v2' else (0,100)
xtitle='v_{2}' if args.xaxis=='v2' else 'q_{2}^{FT0C} percentile (%)'
ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
base=Path(__file__).resolve().parents[3]
out=base/'farm/fit/farm_output/merged_full/oscillation_comparison'
fits=json.loads((out/'fourier_fits_hesse.json').read_text())['fits']
boot=json.loads((out/'FT0C_bootstrap_propagated_errors.json').read_text()) if args.bootstrap_errors else None
original_lines=json.loads((out/'epsilon_vs_v2_inputs.json').read_text())['pol1_fits'] if boot else []
if boot:
 out=out/'FT0C_bootstrap_errors';out.mkdir(parents=True,exist_ok=True)
resolution_file=base/'farm/fit/farm_output/ft0a_resolution_vs_q2.csv'
with resolution_file.open() as f:
 res={int(r['q2_bin']):r for r in csv.DictReader(f) if r['centrality']=='10-30'}
records=[]
linear_results=[]
for v in fits:
 if v['component']!='Rside':continue
 pl=v['plane'];q=v['qn'];r=float(res[q]['res_'+pl]);sr=float(res[q]['err_res'+('A' if pl=='FT0A' else 'C')+'_independent'])
 e=v['epsilon'];se=v['epsilon_error']
 if boot and pl=='FT0C':se=next(b for b in boot['epsilon'] if b['mt']==v['mt'] and b['qn']==q)['error']
 records.append(dict(mt=v['mt'],plane=pl,qn=q,observed=e,observed_error=se,
                     resolution=r,resolution_error=sr,corrected=e/r,
                     corrected_error=math.hypot(se/r,e*sr/r**2)))
for corrected in [False,True]:
 name='epsilon_vs_'+('v2' if args.xaxis=='v2' else 'q2_percentile')+'_'+('resolution_corrected' if corrected else 'observed')+'_vertical'
 value='corrected' if corrected else 'observed'
 c=ROOT.TCanvas(name,name,850,1100);c.Divide(1,2,0,.002)
 keep=[]
 # Same y limits across the two mT panels of a given figure.
 low=min(0.,min(v[value]-v[value+'_error'] for v in records))
 high=max(v[value]+v[value+'_error'] for v in records)
 span=max(.05,high-low)
 for mt in range(2):
  c.cd(mt+1);ROOT.gPad.SetLeftMargin(.16);ROOT.gPad.SetRightMargin(.05);ROOT.gPad.SetBottomMargin(.15);ROOT.gPad.SetTopMargin(.15)
  frame=ROOT.gPad.DrawFrame(xlo,low-.08*span,xhi,high+.42*span)
  frame.SetTitle(';'+xtitle+';'+('#epsilon_{corr} = #epsilon_{obs}/Res{#Psi_{2}}' if corrected else '#epsilon_{obs} = 2R_{side,2}^{2}/R_{side,0}^{2}'))
  if args.xaxis=='v2':frame.GetXaxis().SetNdivisions(505)
  for ax in [frame.GetXaxis(),frame.GetYaxis()]:ax.SetTitleSize(.047);ax.SetLabelSize(.041)
  frame.GetYaxis().SetTitleOffset(1.35)
  leg=ROOT.TLegend(.55,.71,.89,.84);leg.SetBorderSize(0);leg.SetFillStyle(0);leg.SetTextSize(.041)
  for pl,col,marker in [('FT0A',ROOT.kViolet+1,20),('FT0C',ROOT.kOrange+7,25)]:
   g=ROOT.TGraphErrors(3);g.SetName(f'epsilon_{pl}_mt{mt}')
   for q in range(3):
    v=next(x for x in records if x['mt']==mt and x['plane']==pl and x['qn']==q)
    g.SetPoint(q,x[q],v[value]);g.SetPointError(q,ex[q],v[value+'_error'])
   g.SetMarkerStyle(marker);g.SetMarkerSize(1.15);g.SetMarkerColor(col);g.SetLineColor(col);g.SetLineWidth(2);g.Draw('PZ SAME');leg.AddEntry(g,pl+' plane','pe');keep.append(g)
   if args.xaxis=='v2':
    # Center the independent variable for numerical stability; this is pol1.
    line=ROOT.TF1(f'pol1_{name}_{mt}_{pl}','[0]+[1]*(x-0.072)',min(x),max(x))
    line.SetParameters(sum(g.GetPointY(i) for i in range(3))/3.,0.)
    if boot:
     original=next(v for v in original_lines if v['mt']==mt and v['plane']==pl and v['corrected']==corrected)
     slope=original['slope'];slope_error=original['slope_error'] if pl=='FT0A' else next(v for v in boot['slopes'] if v['mt']==mt and v['corrected']==corrected)['error']
     line.SetParameters(original['intercept_at_v2_0p072'],slope)
    else:
     result=g.Fit(line,'SQNR')  # ROOT effective variance includes x errors.
     if int(result)!=0:raise RuntimeError(f'pol1 fit failed: {name} mt{mt} {pl}')
     slope=line.GetParameter(1);slope_error=line.GetParError(1);keep.append(result)
    line.SetLineColor(col);line.SetLineWidth(2);line.Draw('SAME');keep.append(line)
    annotation=ROOT.TLatex();annotation.SetNDC();annotation.SetTextFont(42);annotation.SetTextSize(.040);annotation.SetTextColor(col)
    annotation.DrawLatex(.20,.80 if pl=='FT0A' else .74,f'{pl}: slope = {slope:.2f} #pm {slope_error:.2f}')
    keep.append(annotation)
    linear_results.append(dict(corrected=corrected,mt=mt,plane=pl,slope=slope,slope_error=slope_error,
      intercept_at_v2_0p072=line.GetParameter(0),chi2=None if boot else line.GetChisquare(),ndf=1 if boot else line.GetNDF(),
      note='Original central line retained; only uncertainty updated' if boot else 'Refitted'))
  leg.Draw();keep.append(leg)
  t=ROOT.TLatex();t.SetNDC();t.SetTextFont(42);t.SetTextSize(.041)
  t.DrawLatex(.16,.935,'Pb-Pb 10-30%; m_{T} = '+('0.4-0.7' if mt==0 else '0.7-1.5')+' GeV')
  t.SetTextSize(.027);t.DrawLatex(.16,.875,('FT0C: bootstrap; FT0A: HESSE; '+('resolution corrected' if corrected else 'uncorrected')) if boot else ('HESSE + approximate resolution errors' if corrected else 'HESSE errors; before resolution correction'));keep.append(t)
 c.SaveAs(str(out/(name+'.png')))
 f=ROOT.TFile.Open(str(out/(name+'.root')),'RECREATE');c.Write();f.Close()
(out/'epsilon_resolution_corrected.json').write_text(json.dumps(dict(resolution_source=str(resolution_file),
 assumptions='Same event-weighted resolution for both mT bins; independent epsilon and resolution errors; resolution errors approximate; no finite-phi-bin correction',points=records),indent=2))
for v in records:print(v)
if args.xaxis=='v2':
 (out/'epsilon_vs_v2_inputs.json').write_text(json.dumps(dict(
  centrality='10-30%',q2_percentile_ranges=[[0,30],[30,70],[70,100]],v2=v2,v2_errors=v2err,
  source='User-supplied values; no additional correction applied to v2',
  note='Same v2 coordinates in both mT panels; no epsilon-v2 regression performed',
  points=records,pol1_fits=linear_results,
  fit_note='ROOT effective-variance fit using x/y errors; point correlations unavailable; 3 points and 2 parameters'),indent=2))
 for fit in linear_results:print('POL1',fit)
