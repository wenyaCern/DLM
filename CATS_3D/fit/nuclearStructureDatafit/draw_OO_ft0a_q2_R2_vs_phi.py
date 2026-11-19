#!/usr/bin/env python3
import math,re
from array import array
from pathlib import Path
import ROOT
ROOT.gROOT.SetBatch(True);ROOT.gStyle.SetOptStat(0)
B=Path(__file__).resolve().parents[2]; L=B/'fit/nuclearStructureDatafit/results/oo_ft0a_merged_mT/logs'; O=B/'fit/nuclearStructureDatafit/results/oo_ft0a_merged_mT/summary';O.mkdir(parents=True,exist_ok=True)
phi=[.1,.7,1.3,1.9,2.5,3.1]; names=[(1,'Rout','R_{out}^{2}'),(2,'Rside','R_{side}^{2}'),(3,'Rlong','R_{long}^{2}')]
colors=[ROOT.kRed+1,ROOT.kBlue+1,ROOT.kGreen+2];marks=[20,21,33]; eps=[];keep=[]
def vals(q,i,n):
 r,e=[],[]
 for p in range(6):
  s=(L/f'qn{q}_phi{p}.log').read_text();m=re.findall(rf'^\s*{i}\s+{n}\s+([0-9.eE+-]+)\s+([0-9.eE+-]+)',s,re.M)
  if not m: raise RuntimeError(f'missing {n} q{q} p{p}')
  a,b=map(float,m[-1]);r.append(a*a);e.append(2*a*b)
 return r,e
c=ROOT.TCanvas('cOOFT0A','',1650,540);c.Divide(3,1,.002,.002)
for pad,(idx,n,yl) in enumerate(names,1):
 c.cd(pad);ROOT.gPad.SetTicks();ROOT.gPad.SetBottomMargin(.16);ROOT.gPad.SetLeftMargin(.15)
 graphs=[];allv=[]
 for q in range(3):
  y,ey=vals(q,idx,n);allv += [z for pair in zip(y,ey) for z in (pair[0]-pair[1],pair[0]+pair[1])]
  g=ROOT.TGraphErrors(6,array('d',phi),array('d',y),array('d',[.3]*6),array('d',ey));g.SetMarkerStyle(marks[q]);g.SetMarkerColor(colors[q]);g.SetLineColor(colors[q]);graphs.append((g,y,ey))
 frame=ROOT.gPad.DrawFrame(-.2,min(allv)-.1*(max(allv)-min(allv)),3.4,max(allv)+.18*(max(allv)-min(allv)));frame.SetTitle(f';#Delta#varphi=#varphi_{{pair}}-#Psi_{{2}}^{{FT0A}} (rad);{yl} (fm^{{2}})');keep.append(frame)
 leg=ROOT.TLegend(.67,.68,.93,.9);leg.SetBorderSize(0)
 for q,(g,y,ey) in enumerate(graphs):
  g.Draw('P SAME');f=ROOT.TF1(f'f_{n}_{q}','[0]+[1]*cos(2*x)',-.2,3.4);f.SetParameters(sum(y)/6,.1);f.SetLineColor(colors[q]);f.SetLineWidth(3);g.Fit(f,'Q0');f.Draw('SAME');leg.AddEntry(g,f'q_{{2}} bin {q}','pe');keep += [g,f]
  if idx==2:
   A,Bb,eA,eB=f.GetParameter(0),f.GetParameter(1),f.GetParError(0),f.GetParError(1);eps.append((q,Bb/A,math.sqrt((eB/A)**2+(Bb*eA/A**2)**2)))
 leg.Draw();keep.append(leg)
c.SaveAs(str(O/'OO_FT0A_cent10to30_mt0p244to0p714_R2_vs_phi.png'));c.SaveAs(str(O/'OO_FT0A_cent10to30_mt0p244to0p714_R2_vs_phi.pdf'))
ce=ROOT.TCanvas('cOOeps','',700,600);x=array('d',[q for q,_,_ in eps]);y=array('d',[v for _,v,_ in eps]);ey=array('d',[e for *_,e in eps]);g=ROOT.TGraphErrors(3,x,y,array('d',[0]*3),ey);g.SetTitle(';q_{2} bin;#epsilon_{final}=2R_{side,2}^{2}/R_{side,0}^{2}');g.SetMarkerStyle(20);g.SetMarkerSize(1.4);g.SetLineWidth(2);g.Draw('APL');ce.SaveAs(str(O/'OO_FT0A_epsilon_vs_q2.png'));ce.SaveAs(str(O/'OO_FT0A_epsilon_vs_q2.pdf'))
for q,v,e in eps: print(f'q2 {q}: epsilon={v:.6f} +/- {e:.6f}')
