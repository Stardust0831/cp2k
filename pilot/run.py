from pathlib import Path
from datetime import datetime,timezone
import hashlib,json,os,subprocess,time
root=Path(__file__).resolve().parent.parent
d=root/'diagnostics-full-spectrum'
assert (root/'source/REVISION').read_text().strip()=='5d4c07bd2d03f4106c5e72172c7a20d5e697dcdd'
art=json.loads((root/'artifacts/baseline/hashes.json').read_text())
for name,h in art.items():
 if name.startswith('build/'):
  assert hashlib.sha256((root/name).read_bytes()).hexdigest()==h,name
result=dict(job=os.environ['SLURM_JOB_ID'],host=os.uname().nodename,source='5d4c07bd2d03f4106c5e72172c7a20d5e697dcdd',artifacts=art,cases={})
mpi=['mpirun','-np','4','--map-by','ppr:4:node:pe=4','--bind-to','core','--report-bindings']
for name in ('real-gth','real-sr','complex-sr'):
 print('START',name,flush=True)
 w=d/name;t=time.monotonic()
 with (w/'cp2k.out').open('w') as f:
  p=subprocess.run(mpi+['python3',str(root/'rank.py'),'-i','input.inp'],cwd=w,stdout=f,stderr=subprocess.STDOUT,timeout=180)
 result['cases'][name]=dict(exit_code=p.returncode,seconds=time.monotonic()-t,input_sha256=hashlib.sha256((w/'input.inp').read_bytes()).hexdigest())
 print('END',name,result['cases'][name],flush=True)
result['finished_at']=datetime.now(timezone.utc).isoformat()
(d/'result.json').write_text(json.dumps(result,indent=2)+'\n')
assert all(result['cases'][name]['exit_code']==0 for name in ('real-sr','complex-sr'))
