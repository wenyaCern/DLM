"""Overlay nominal FT0A/FT0C R^2 vs phi with propagated HESSE errors.
Weighted cosine fits assume independent phi points; include A/B covariance in epsilon.
Phi coordinates are raw merged-axis centers; no resolution correction applied.
"""
from pathlib import Path
import math
import json
import argparse
import ROOT
parser=argparse.ArgumentParser()
parser.add_argument('--bootstrap-errors',action='store_true')
args=parser.parse_args()
ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
BASE=Path(__file__).resolve().parents[3]
OUT=BASE/'farm/fit/farm_output/merged_full/oscillation_comparison'
reference=OUT
boot=json.loads((reference/'FT0C_bootstrap_propagated_errors.json').read_text()) if args.bootstrap_errors else None
original_fits=json.loads((reference/'fourier_fits_hesse.json').read_text())['fits'] if boot else []
if boot:OUT=reference/'FT0C_bootstrap_errors'
OUT.mkdir(parents=True,exist_ok=True)
phi=[.1,.7,1.3,1.9,2.5,3.1]
rows=[]; data={}
for mt in range(2):
 for plane in ['FT0A','FT0C']:
  for q in range(3):
   for p in range(6):
    path=BASE/f'farm/fit/farm_output/merged_full/{plane}/cent1/mtMerged{mt}/qn{q}/fit_phi{p}.root'
    f=ROOT.TFile.Open(str(path))
    if not f or f.IsZombie(): raise RuntimeError(str(path))
    t=f.Get('fit_summary')
    if not t or not t.GetEntries():raise RuntimeError(f'Missing summary {path}')
    t.GetEntry(0)
    record=dict(mt=mt,plane=plane,qn=q,phi_bin=p,phi_center_rad=phi[p])
    for r in ['Rout','Rside','Rlong','Chi2Ndf','FitStatus']:
     record[r]=float(t.GetLeaf(r).GetValue())
    for r in ['Rout','Rside','Rlong']:
     e=f.Get(r+'_error')
     record[r+'Err']=float(e.GetVal()) if e else float('nan')
    for key in ['MIGRAD_status','baseline_MIGRAD_status','HESSE_status']:
     obj=f.Get(key);record[key]=int(obj.GetVal()) if obj else None
    if record['MIGRAD_status']!=0 or record['baseline_MIGRAD_status']!=0:
     raise RuntimeError(f'Inspect fit convergence before plotting: {record}')
    if any(not math.isfinite(record[r]) or record[r]<=0 or not math.isfinite(record[r+'Err']) for r in ['Rout','Rside','Rlong']):
     raise RuntimeError(f'Invalid radius {path}')
    data[mt,plane,q,p]=record;rows.append(record);f.Close()

def cosine_fit(y,e):
 if any(v<=0 for v in e):raise RuntimeError('Nonpositive HESSE error')
 z=[math.cos(2*x) for x in phi];w=[1/v**2 for v in e]
 s=sum(w);sz=sum(a*b for a,b in zip(w,z));szz=sum(a*b*b for a,b in zip(w,z))
 sy=sum(a*b for a,b in zip(w,y));szy=sum(a*b*c for a,b,c in zip(w,z,y));det=s*szz-sz*sz
 A=(szz*sy-sz*szy)/det;B=(s*szy-sz*sy)/det
 va=szz/det;vb=s/det;cov=-sz/det
 eps=B/A;variance=vb/A**2+B**2*va/A**4-2*B*cov/A**3
 chi2=sum(((v-A-B*x)/err)**2 for v,x,err in zip(y,z,e))
 return dict(A=A,B=B,R2_0=A,R2_2=B/2,A_error=math.sqrt(va),B_error=math.sqrt(vb),covAB=cov,
             epsilon=eps,epsilon_error=math.sqrt(max(0.,variance)),chi2=chi2,ndf=4)

fits=[]

for mt in range(2):
 for squared in [True]:
  name=('R2' if squared else 'R')+f'_vs_phi_FT0A_FT0C_mtMerged{mt}'
  c=ROOT.TCanvas(name,name,1500,1080)
  header=ROOT.TPad('header_'+name,'',0,.91,1,1);header.Draw();header.cd()
  title=ROOT.TLatex();title.SetNDC();title.SetTextFont(42);title.SetTextSize(.27)
  mtlabel='0.4-0.7' if mt==0 else '0.7-1.5'
  title.DrawLatex(.06,.62,f'Pb-Pb 10-30%; m_{{T}} = {mtlabel} GeV; FT0C q_{{2}} selection')
  title.SetTextSize(.21);title.DrawLatex(.06,.20,'FT0C: bootstrap; FT0A: HESSE; original fits retained; no resolution correction' if boot else 'HESSE errors; weighted A+B cos(2#Delta#varphi) fits; no event-plane resolution correction')
  c.cd();body=ROOT.TPad('body_'+name,'',0,0,1,.91);body.Draw();body.Divide(3,3,.001,.001)
  keep=[header,body,title]
  for rindex,r in enumerate(['Rout','Rside','Rlong']):
   for q in range(3):
    body.cd(3*rindex+q+1);ROOT.gPad.SetLeftMargin(.16);ROOT.gPad.SetBottomMargin(.19);ROOT.gPad.SetTopMargin(.08)
    ys={pl:[data[mt,pl,q,p][r]**(2 if squared else 1) for p in range(6)] for pl in ['FT0A','FT0C']}
    es={pl:[(2*data[mt,pl,q,p][r]*data[mt,pl,q,p][r+'Err'] if squared else data[mt,pl,q,p][r+'Err']) for p in range(6)] for pl in ['FT0A','FT0C']}
    curves={pl:cosine_fit(y,es[pl]) for pl,y in ys.items()}
    if boot:
     es['FT0C']=[next(v for v in boot['points'] if v['mt']==mt and v['qn']==q and v['phi']==p)['R2_errors'][r] for p in range(6)]
     curves={pl:next(v for v in original_fits if v['mt']==mt and v['qn']==q and v['plane']==pl and v['component']==r) for pl in ys}
    for pl,result in curves.items():fits.append({**result,'mt':mt,'plane':pl,'qn':q,'component':r})
    allvals=sum(ys.values(),[])+[v+s for pl in es for v,s in zip(ys[pl],es[pl])]+[v-s for pl in es for v,s in zip(ys[pl],es[pl])]
    if squared:
     allvals += [v['A']+sgn*abs(v['B']) for v in curves.values() for sgn in [-1,1]]
    lo=min(allvals);hi=max(allvals);pad=max(.12 if squared else .025,(hi-lo)*.18)
    frame=ROOT.gPad.DrawFrame(-.2,lo-pad,3.4,hi+pad*2.5)
    label=r[1:];unit='fm^{2}' if squared else 'fm';power='^{2}' if squared else ''
    frame.SetTitle(f';#Delta#varphi (rad);R_{{{label}}}{power} ({unit})')
    for ax in [frame.GetXaxis(),frame.GetYaxis()]:ax.SetTitleSize(.064);ax.SetLabelSize(.05)
    frame.GetYaxis().SetTitleOffset(1.12)
    legend=ROOT.TLegend(.22,.77,.88,.90);legend.SetNColumns(2);legend.SetBorderSize(0);legend.SetFillStyle(0);legend.SetTextSize(.049)
    for pl,col,marker in [('FT0A',ROOT.kViolet+1,20),('FT0C',ROOT.kOrange+7,25)]:
     g=ROOT.TGraphErrors(6)
     for p,y in enumerate(ys[pl]):g.SetPoint(p,phi[p],y);g.SetPointError(p,0,es[pl][p])
     g.SetMarkerStyle(marker);g.SetMarkerSize(.95);g.SetMarkerColor(col);g.SetLineColor(col);g.Draw('P SAME');keep.append(g)
     legend.AddEntry(g,pl+' plane','p')
     if squared:
      a,b=curves[pl]['A'],curves[pl]['B'];line=ROOT.TF1(f'cos_{mt}_{r}_{q}_{pl}','[0]+[1]*cos(2*x)',-.2,3.4);line.SetParameters(a,b);line.SetLineColor(col);line.SetLineWidth(2);line.SetLineStyle(2);line.Draw('SAME');keep.append(line)
    legend.Draw();keep.append(legend)
    text=ROOT.TLatex();text.SetNDC();text.SetTextSize(.052);text.DrawLatex(.20,.68,f'q_{{2}} bin {q}');keep.append(text)
  c.SaveAs(str(OUT/(name+'.png')))
  f=ROOT.TFile.Open(str(OUT/(name+'.root')),'RECREATE');c.Write();f.Close()
(OUT/'nominal_points.json').write_text(json.dumps(rows,indent=2))
print('Validated 72 fits. Wrote R and R2 overlays for both merged mT bins:',OUT)
if boot:raise SystemExit(0)  # Eccentricity is rendered by draw_epsilon_vertical.py.

(OUT/'fourier_fits_hesse.json').write_text(json.dumps(dict(
 convention='R^2=A+B*cos(2phi); epsilon=B/A from side only',
 assumptions='Diagonal HESSE point covariance; raw phi-bin centers; no finite-bin or EP resolution correction',
 fits=fits),indent=2))
c=ROOT.TCanvas('eccentricity','Eccentricity vs q2 percentile',1400,650);c.Divide(2,1)
keep=[]
for mt in range(2):
 c.cd(mt+1);ROOT.gPad.SetLeftMargin(.15);ROOT.gPad.SetBottomMargin(.15);ROOT.gPad.SetTopMargin(.14)
 selected=[v for v in fits if v['mt']==mt and v['component']=='Rside']
 lo=min(0.,min(v['epsilon']-v['epsilon_error'] for v in selected))
 hi=max(v['epsilon']+v['epsilon_error'] for v in selected);span=max(.05,hi-lo)
 frame=ROOT.gPad.DrawFrame(0,lo-.1*span,100,hi+.5*span)
 frame.SetTitle(';q_{2}^{FT0C} percentile (%);#epsilon_{obs} = 2R_{side,2}^{2}/R_{side,0}^{2}')
 for ax in [frame.GetXaxis(),frame.GetYaxis()]:ax.SetTitleSize(.042);ax.SetLabelSize(.036)
 frame.GetYaxis().SetTitleOffset(1.6)
 legend=ROOT.TLegend(.50,.72,.89,.85);legend.SetBorderSize(0);legend.SetFillStyle(0)
 for pl,col,marker in [('FT0A',ROOT.kViolet+1,20),('FT0C',ROOT.kOrange+7,25)]:
  g=ROOT.TGraphErrors(3)
  for q in range(3):
   v=next(v for v in selected if v['plane']==pl and v['qn']==q)
   g.SetPoint(q,[15,50,85][q],v['epsilon']);g.SetPointError(q,[15,20,15][q],v['epsilon_error'])
   print(f'mT{mt} {pl} q{q}: epsilon={v["epsilon"]:.6f} +/- {v["epsilon_error"]:.6f}; Fourier chi2/4={v["chi2"]:.3f}/4')
  g.SetMarkerStyle(marker);g.SetMarkerSize(1.2);g.SetMarkerColor(col);g.SetLineColor(col);g.SetLineWidth(2);g.Draw('PZ SAME');legend.AddEntry(g,pl+' plane','pe');keep.append(g)
 legend.Draw();keep.append(legend)
 t=ROOT.TLatex();t.SetNDC();t.SetTextFont(42);t.SetTextSize(.034)
 t.DrawLatex(.15,.92,'Pb-Pb 10-30%; m_{T} = '+('0.4-0.7' if mt==0 else '0.7-1.5')+' GeV')
 t.SetTextSize(.027);t.DrawLatex(.15,.87,'HESSE propagation; no EP-resolution correction');keep.append(t)
c.SaveAs(str(OUT/'epsilon_vs_q2_percentile_hesse.png'));c.SaveAs(str(OUT/'epsilon_vs_q2_percentile_hesse.pdf'))
f=ROOT.TFile.Open(str(OUT/'epsilon_vs_q2_percentile_hesse.root'),'RECREATE');c.Write();f.Close()
