"""Compare only overlapping intervals, recording incomplete refinement explicitly."""
from pathlib import Path
import csv, json, tomllib
import numpy as np

root=Path(__file__).resolve().parents[1]
def load(name):
    p=root/'results'/name
    with (p/'history.csv').open() as f:
        r=csv.DictReader(f);keys=r.fieldnames or ['time']
        rows=[q for q in r if all(q.get(k) is not None for k in keys)]
    return ({k:np.array([float(q[k]) for q in rows]) for k in keys},
            tomllib.loads((p/'run.toml').read_text()))

pairs=[('face_norm_choice','bubble-face-norm','bubble-proposed-initial'),('time_7p5_15','bubble-proposed-factor7p5','bubble-proposed-factor15-fixed'),('standard_time_half','bubble-standard-half','bubble-standard-mg'),('grid_128_256','bubble-proposed-128','bubble-proposed-initial'),
       ('grid_256_512','bubble-proposed-initial','bubble-proposed-512-fixed'),
       ('time_15_30','bubble-proposed-factor15-fixed','bubble-proposed-initial'),
       ('roundoff_fix','bubble-proposed-factor15','bubble-proposed-factor15-fixed'),
       ('backend_jacobi_multigrid','bubble-standard-initial','bubble-standard-mg')]
out={}
for label,a,b in pairs:
    x,mx=load(a);y,my=load(b)
    if min(len(x['time']),len(y['time']))<2:continue
    end=min(x['time'][-1],y['time'][-1],.07)
    t=np.linspace(0,end,1001)
    diff=np.interp(t,x['time'],x['rise_velocity'])-np.interp(t,y['time'],y['rise_velocity'])
    out[label]={'runs':[a,b],'end_s':float(end), 'full_interval':bool(end>=.07-1e-12),
                'rms_difference_m_s':float(np.sqrt(np.sum((diff[:-1]**2+diff[1:]**2)*np.diff(t)/2)/end)),
                'max_difference_m_s':float(max(abs(diff))),
                'source_hashes':[mx['source_sha256'],my['source_sha256']],
                'note':'No Richardson/GCI estimate; epsilon=dx changes with the grid.'}
p=root/'evidence/local-20260907/sensitivity.json'
p.write_text(json.dumps(out,indent=2)+'\n');print(p.read_text())
