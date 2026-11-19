"""Propagate paired replicas through the original fixed-weight estimators.
Nominal Fourier/pol1 central values are not refitted. Phi/q bins must share draws.
"""
from pathlib import Path
import json,math
from statistics import stdev
import ROOT
ROOT.gROOT.SetBatch(True)
base=Path(__file__).resolve().parents[3]
old=base/'farm/fit/farm_output/merged_full/oscillation_comparison'
points=json.loads((old/'nominal_points.json').read_text())
fits=json.loads((old/'fourier_fits_hesse.json').read_text())['fits']
saved=json.loads((old/'epsilon_vs_v2_inputs.json').read_text())
z=[math.cos(2*x) for x in [.1,.7,1.3,1.9,2.5,3.1]]
def line(x,y,w):
 s=sum(w);sx=sum(a*b for a,b in zip(x,w));sxx=sum(a*a*b for a,b in zip(x,w))
 sy=sum(a*b for a,b in zip(y,w));sxy=sum(a*b*c for a,b,c in zip(x,y,w));det=s*sxx-sx*sx
 return (sxx*sy-sx*sxy)/det,(s*sxy-sx*sy)/det
rep={};draws={};result=dict(points=[],epsilon=[],slopes=[],method='Paired bootstrap; original HESSE estimator weights fixed; nominal centers and lines unchanged')
for mt in range(2):
 for q in range(3):
  for p in range(6):
   values=[]
   for b in range(100):
    path=base/f'farm/fit/farm_output/merged_bootstrap/FT0C/cent1/mtMerged{mt}/qn{q}/phi{p}/fit_{b:03d}.root'
    f=ROOT.TFile.Open(str(path))
    if not f or f.IsZombie():raise RuntimeError(str(path))
    for key in ['MIGRAD_status','baseline_MIGRAD_status']:
     if not f.Get(key) or f.Get(key).GetVal()!=0:raise RuntimeError(f'{key}: {path}')
    d=f.Get('bootstrap_draw_indices').GetTitle()
    if b in draws and draws[b]!=d:raise RuntimeError(f'Unpaired replica {path}')
    draws[b]=d;t=f.Get('fit_summary');t.GetEntry(0)
    values.append({r:float(t.GetLeaf(r).GetValue()) for r in ['Rout','Rside','Rlong']});f.Close()
   rep[mt,q,p]=values
   result['points'].append(dict(mt=mt,qn=q,phi=p,R2_errors={r:stdev(v[r]**2 for v in values) for r in ['Rout','Rside','Rlong']}))
   print(f'Read mt{mt} q{q} phi{p}',flush=True)
eps={}
for mt in range(2):
 for q in range(3):
  nominal=[next(v for v in points if v['mt']==mt and v['plane']=='FT0C' and v['qn']==q and v['phi_bin']==p) for p in range(6)]
  w=[1/(2*v['Rside']*v['RsideErr'])**2 for v in nominal]
  e=[]
  for b in range(100):
   A,B=line(z,[rep[mt,q,p][b]['Rside']**2 for p in range(6)],w);e.append(B/A)
  eps[mt,q]=e;result['epsilon'].append(dict(mt=mt,qn=q,error=stdev(e)))
 for corrected in [False,True]:
  previous=next(v for v in saved['pol1_fits'] if v['mt']==mt and v['plane']=='FT0C' and v['corrected']==corrected)
  records=[next(v for v in saved['points'] if v['mt']==mt and v['plane']=='FT0C' and v['qn']==q) for q in range(3)]
  key='corrected' if corrected else 'observed'
  w=[1/(v[key+'_error']**2+previous['slope']**2*saved['v2_errors'][q]**2) for q,v in enumerate(records)]
  x=[v-.072 for v in saved['v2']];slopes=[]
  for b in range(100):
   y=[eps[mt,q][b]/(records[q]['resolution'] if corrected else 1.) for q in range(3)]
   slopes.append(line(x,y,w)[1])
  # Independent external x and resolution errors, linearized with fixed weights.
  s=sum(w);sx=sum(a*b for a,b in zip(x,w));sxx=sum(a*a*b for a,b in zip(x,w));det=s*sxx-sx*sx
  coeff=[ww*(s*xx-sx)/det for xx,ww in zip(x,w)]
  extra=sum((coeff[q]*previous['slope']*saved['v2_errors'][q])**2 for q in range(3))
  if corrected:extra+=sum((coeff[q]*v['observed']*v['resolution_error']/v['resolution']**2)**2 for q,v in enumerate(records))
  result['slopes'].append(dict(mt=mt,corrected=corrected,error=math.sqrt(stdev(slopes)**2+extra)))
(old/'FT0C_bootstrap_propagated_errors.json').write_text(json.dumps(result,indent=2))
print('Saved paired bootstrap errors for 36 bins, six eccentricities and four slopes.')
