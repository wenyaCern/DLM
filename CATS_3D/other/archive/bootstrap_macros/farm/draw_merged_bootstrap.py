"""Local PyROOT plotting: nominal centers and pointwise bootstrap standard deviations.
No HESSE errors or pre-existing analytic bands are used. All requested replicas
must exist and have successful MIGRAD; failures stop instead of silently dropping.
"""
import argparse
import json
import math
from pathlib import Path
from statistics import stdev, mean
import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
p=argparse.ArgumentParser()
p.add_argument('nominal')
p.add_argument('replica_directory')
p.add_argument('output_directory')
p.add_argument('--count',type=int,default=100)
p.add_argument('--projection',choices=['unweighted','mixed'],default='unweighted')
p.add_argument('--label',default='')
a=p.parse_args()
if a.count<2: raise RuntimeError('Need at least two replicas')
out=Path(a.output_directory)
if out.exists(): raise RuntimeError('Use a new output directory; no overwrite')
keys=['hDataCk3D','hCATS3D_background','hCATS3D_best']
suffix='_proj_' if a.projection=='unweighted' else '_mixedWeighted_proj_'

def read(path,replica=False):
    f=ROOT.TFile.Open(str(path))
    if not f or f.IsZombie(): raise RuntimeError(f'Missing {path}')
    status=f.Get('MIGRAD_status')
    if not status or status.GetVal()!=0: raise RuntimeError(f'Inspect MIGRAD failure/missing status: {path}')
    baseline=f.Get('baseline_MIGRAD_status')
    if not baseline or baseline.GetVal()!=0: raise RuntimeError(f'Inspect baseline prefit status: {path}')
    t=f.Get('fit_summary'); t.GetEntry(0)
    radii={r:float(t.GetLeaf(r).GetValue()) for r in ['Rout','Rside','Rlong']}
    if any(not math.isfinite(v) or v<=0 for v in radii.values()):raise RuntimeError(f'Invalid radius: {path}')
    if int(t.GetLeaf('NFitBins').GetValue())!=1000:raise RuntimeError(f'Unexpected fit range/bin count: {path}')
    hs={}
    for d in ['out','side','long']:
        for k in keys:
            name=k+suffix+d
            h=f.Get(name)
            if not h: raise RuntimeError(f'Missing {name}: {path}')
            h=h.Clone();h.SetDirectory(0);hs[name]=h
    draws=f.Get('bootstrap_draw_indices')
    if replica and not draws: raise RuntimeError(f'Missing draw provenance: {path}')
    cut=f.Get('slice_cut_abs')
    if not cut: raise RuntimeError(f'Missing slice metadata: {path}')
    result=(hs,radii,float(cut.GetVal()),int(f.Get('HESSE_status').GetVal()))
    f.Close();return result

nom,nr,cut,nh=read(a.nominal)
replicas=[read(Path(a.replica_directory)/f'fit_{i:03d}.root',True) for i in range(a.count)]
for hs,rs,c,h in replicas:
    if c!=cut: raise RuntimeError('Different slice cuts')
    for name,ref in nom.items():
        x=hs[name]
        if x.GetNbinsX()!=ref.GetNbinsX() or any(abs(x.GetXaxis().GetBinLowEdge(b)-ref.GetXaxis().GetBinLowEdge(b))>1e-9 for b in range(1,ref.GetNbinsX()+2)):
            raise RuntimeError('Different projection axes')
        if any(not math.isfinite(x.GetBinContent(b)) for b in range(1,x.GetNbinsX()+1)):
            raise RuntimeError('Nonfinite replica projection')
out.mkdir(parents=True)
summary=dict(count=a.count,projection=a.projection,slice_cut_abs=cut,
             nominal=a.nominal,hesse_nonzero=sum(r[3]!=0 for r in replicas),
             errors='pointwise sample standard deviation (N-1); period-block bootstrap',radii={},residuals={})
for r in nr:
    values=[x[1][r] for x in replicas]
    summary['radii'][r]=dict(nominal=nr[r],bootstrap_mean=mean(values),bootstrap_sigma=stdev(values))
root=ROOT.TFile.Open(str(out/'bootstrap_bands.root'),'CREATE')
for d in ['out','side','long']:
    graphs=[]
    for k in keys:
        name=k+suffix+d;h=nom[name];g=ROOT.TGraphErrors(h.GetNbinsX());g.SetName(name+'_bootstrap')
        for b in range(1,h.GetNbinsX()+1):
            values=[x[0][name].GetBinContent(b) for x in replicas]
            g.SetPoint(b-1,h.GetBinCenter(b),h.GetBinContent(b))
            g.SetPointError(b-1,h.GetBinWidth(b)/2,stdev(values))
        graphs.append(g);g.Write()
    data,bkg,fit=graphs
    c=ROOT.TCanvas('c_'+d,'Bootstrap bands',950,900)
    top=ROOT.TPad('top_'+d,'',0,.29,1,1);bottom=ROOT.TPad('bottom_'+d,'',0,0,1,.29)
    top.SetLeftMargin(.14);top.SetBottomMargin(.025);top.SetTopMargin(.07)
    bottom.SetLeftMargin(.14);bottom.SetTopMargin(.02);bottom.SetBottomMargin(.29)
    top.Draw();bottom.Draw();top.cd()
    frame=nom[keys[0]+suffix+d].Clone('frame_'+d);frame.Reset()
    vals=[(g.GetPointY(i)-g.GetErrorY(i),g.GetPointY(i)+g.GetErrorY(i)) for g in graphs for i in range(g.GetN())]
    lo=min(x[0] for x in vals);hi=max(x[1] for x in vals);pad=max(.05,(hi-lo)*.12)
    frame.SetMinimum(lo-pad);frame.SetMaximum(hi+max(pad,(hi-lo)*.32))
    frame.SetTitle('; ;C(q)');frame.GetXaxis().SetLabelSize(0);frame.GetYaxis().SetTitleSize(.05);frame.GetYaxis().SetTitleOffset(1.2);frame.Draw('AXIS')
    bkg.SetFillColorAlpha(ROOT.kOrange+1,.3);bkg.SetLineColor(ROOT.kOrange+7);bkg.Draw('2 SAME')
    fit.SetFillColorAlpha(ROOT.kBlue,.22);fit.SetLineColor(ROOT.kBlue);fit.Draw('2 SAME')
    bkgline=nom[keys[1]+suffix+d];bkgline.SetLineColor(ROOT.kOrange+7);bkgline.SetLineWidth(2);bkgline.Draw('HIST SAME')
    fitline=nom[keys[2]+suffix+d];fitline.SetLineColor(ROOT.kBlue);fitline.SetLineWidth(2);fitline.Draw('HIST SAME')
    data.SetMarkerStyle(20);data.SetMarkerSize(.8);data.Draw('PZ SAME')
    leg=ROOT.TLegend(.51,.73,.88,.89);leg.SetBorderSize(0);leg.SetFillStyle(0)
    leg.AddEntry(data,'Data: bootstrap errors','pe');leg.AddEntry(fit,'Nominal fit + bootstrap 1#sigma','lf');leg.AddEntry(bkg,'Baseline + bootstrap 1#sigma','lf');leg.Draw()
    label=ROOT.TLatex();label.SetNDC();label.SetTextSize(.026)
    label.DrawLatex(.15,.95,a.label or (f'{a.count} bootstrap replicas; unweighted slice' if a.projection=='unweighted' else f'{a.count} bootstrap replicas; mixed-weighted slice'))
    bottom.cd()
    residual=ROOT.TGraphErrors();residual.SetName('paired_bootstrap_standardized_residual_'+d)
    residual_rows=[]
    for b in range(1,nom[keys[0]+suffix+d].GetNbinsX()+1):
        differences=[x[0][keys[0]+suffix+d].GetBinContent(b)-x[0][keys[2]+suffix+d].GetBinContent(b) for x in replicas]
        sigma=stdev(differences)
        delta=nom[keys[0]+suffix+d].GetBinContent(b)-nom[keys[2]+suffix+d].GetBinContent(b)
        if sigma<=0:raise RuntimeError(f'Zero paired-residual variance {d} bin {b}')
        h=nom[keys[0]+suffix+d];i=residual.GetN()
        residual.SetPoint(i,h.GetBinCenter(b),delta/sigma);residual.SetPointError(i,h.GetBinWidth(b)/2,0.)
        residual_rows.append(dict(q=h.GetBinCenter(b),delta=delta,sigma=sigma,z=delta/sigma,bootstrap_mean_residual=mean(differences)))
    summary['residuals'][d]=residual_rows
    limit=max(3.5,max(abs(v['z']) for v in residual_rows)*1.15)
    lower=bottom.DrawFrame(frame.GetXaxis().GetXmin(),-limit,frame.GetXaxis().GetXmax(),limit)
    lower.SetTitle(f';q_{{{d}}} (GeV/c);(data-fit)/#sigma_{{boot}}')
    for ax in [lower.GetXaxis(),lower.GetYaxis()]:ax.SetTitleSize(.10);ax.SetLabelSize(.085)
    lower.GetYaxis().SetTitleOffset(.60);lower.GetYaxis().SetNdivisions(505)
    lines=[]
    for y in [-3,-2,-1,0,1,2,3]:
        line=ROOT.TLine(frame.GetXaxis().GetXmin(),y,frame.GetXaxis().GetXmax(),y)
        line.SetLineColor(ROOT.kGray+1);line.SetLineStyle(1 if y==0 else 3);line.Draw();lines.append(line)
    residual.SetMarkerStyle(20);residual.SetMarkerSize(.8);residual.Draw('PZ SAME');root.cd();residual.Write()
    c.cd()
    c.SaveAs(str(out/f'slice_{d}.png'));c.Write()
ROOT.TNamed('bootstrap_summary',json.dumps(summary)).Write();root.Close()
(out/'bootstrap_summary.json').write_text(json.dumps(summary,indent=2))
print(json.dumps(summary,indent=2))
