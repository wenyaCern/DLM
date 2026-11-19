"""Render all 36 FT0C bins with 100 replicas each; preserve failures explicitly."""
from pathlib import Path
import subprocess
import sys
import json
from concurrent.futures import ThreadPoolExecutor,as_completed
base=Path(__file__).resolve().parents[3]
out=base/'other/figures/bootstrap_FT0C_all_merged'
out.mkdir(parents=True,exist_ok=True)
def run(mt,q,p):
 dest=out/f'mtMerged{mt}/qn{q}/phi{p}'
 nominal=base/f'farm/fit/farm_output/merged_full/FT0C/cent1/mtMerged{mt}/qn{q}/fit_phi{p}.root'
 source=base/f'farm/fit/farm_output/merged_bootstrap/FT0C/cent1/mtMerged{mt}/qn{q}/phi{p}'
 label=f'FT0C; 10-30%; m_{{T}} '+('0.4-0.7' if mt==0 else '0.7-1.5')+f' GeV; q_{{2}} bin {q}; #varphi bin {p}; 100 replicas'
 log=out/f'mt{mt}_qn{q}_phi{p}.log'
 if dest.exists():raise RuntimeError(f'Output exists: {dest}')
 with log.open('w') as f:
  r=subprocess.run([sys.executable,str(base/'farm/draw_merged_bootstrap.py'),str(nominal),str(source),str(dest),'--label',label],stdout=f,stderr=subprocess.STDOUT)
 if r.returncode:return dict(mt=mt,qn=q,phi=p,error=str(log))
 s=json.loads((dest/'bootstrap_summary.json').read_text())
 return dict(mt=mt,qn=q,phi=p,count=s['count'],hesse_nonzero=s['hesse_nonzero'],radii=s['radii'],directory=str(dest))
results=[]
with ThreadPoolExecutor(max_workers=2) as executor:
 tasks=[executor.submit(run,mt,q,p) for mt in range(2) for q in range(3) for p in range(6)]
 for t in as_completed(tasks):
  v=t.result();results.append(v);print(f'{len(results)}/36: mt{v["mt"]} qn{v["qn"]} phi{v["phi"]}: '+('FAILED' if 'error' in v else 'OK'),flush=True)
results.sort(key=lambda v:(v['mt'],v['qn'],v['phi']))
(out/'R_uncertainty_summary.json').write_text(json.dumps(results,indent=2))
lines=['# FT0C merged bootstrap radius uncertainties','',
       'Centrality 10–30%; 100 replicas/bin. Values are nominal R ± bootstrap sample standard deviation, in fm.','',
       '| mT (GeV) | q2 bin | phi bin | Rout (fm) | Rside (fm) | Rlong (fm) |',
       '|---|---:|---:|---:|---:|---:|']
for v in results:
 if 'error' in v:
  lines.append(f'| {v["mt"]} | {v["qn"]} | {v["phi"]} | FAILED | See log | |');continue
 vals=[f'{v["radii"][r]["nominal"]:.4f} ± {v["radii"][r]["bootstrap_sigma"]:.4f}' for r in ['Rout','Rside','Rlong']]
 lines.append('| '+('0.4–0.7' if v['mt']==0 else '0.7–1.5')+f' | {v["qn"]} | {v["phi"]} | '+' | '.join(vals)+' |')
(out/'R_uncertainty_summary.md').write_text('\n'.join(lines)+'\n')
print(f'Complete: {sum("error" not in v for v in results)}/36 bins; {out}',flush=True)
