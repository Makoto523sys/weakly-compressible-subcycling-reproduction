"""Independent final shape review; requires complete runs at the same physical time."""
from pathlib import Path
import csv,json,tomllib
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'evidence/local-20260907'

def field(run):
    p=ROOT/'results'/run
    meta=tomllib.loads((p/'run.toml').read_text())
    if meta['status']!='completed_unvalidated':
        raise RuntimeError(f'{run} is incomplete; no final shape comparison')
    index=p/'snapshots.csv'
    if not index.exists():index=p/'snapshots-derived.csv'
    with index.open() as f:last=list(csv.DictReader(f))[-1]
    assert abs(float(last['time'])-.07)<1e-12
    c=meta['config'];nx,ny=c['nx'],c['ny']
    text=(p/last['file']).read_text()
    phi=np.fromstring(text.split('LOOKUP_TABLE default\n')[1].split('SCALARS')[0],sep=' ').reshape((nx,ny),order='F')
    return c,phi,last

names=['bubble-proposed-128','bubble-proposed-initial','bubble-proposed-512-fixed']
fig,ax=plt.subplots(figsize=(6,8));stats={};contours={};arrays={}
for name in names:
    c,phi,last=field(name);nx,ny=c['nx'],c['ny'];dx=c['lx']/nx;dy=c['ly']/ny
    x=(np.arange(nx)+.5)*dx;y=(np.arange(ny)+.5)*dy
    cs=ax.contour(x*1000,y*1000,phi.T,levels=[.5],colors=[f'C{len(stats)}'])
    paths=[p for p in cs.allsegs[0] if len(p)>2]
    pts=np.concatenate(paths)/1000;contours[nx]=pts;arrays[nx]=phi
    gas=1-phi;area=gas.sum()*dx*dy
    stats[name]={'snapshot':last,'gas_area_m2':float(area),
      'centroid_x_m':float(np.sum(gas*x[:,None])/gas.sum()),
      'centroid_y_m':float(np.sum(gas*y[None,:])/gas.sum()),
      'mirror_phase_linf':float(np.max(np.abs(phi-phi[::-1,:]))),
      'contour_components':len(paths),'phi_min':float(phi.min()),'phi_max':float(phi.max())}
    for k,p in enumerate(paths):
        np.savetxt(OUT/f'contour-{nx}-{k}.csv',p/1000,delimiter=',',header='x_m,y_m',comments='')
    ax.plot([],[],color=f'C{len(stats)-1}',label=f'{nx} x {ny}')

def directed(a,b):
    return max(float(np.sqrt(np.min(np.sum((q[:,None,:]-b[None,:,:])**2,axis=2),axis=1)).max()) for q in np.array_split(a,max(1,len(a)//128)))
for coarse,fine in [(128,256),(256,512)]:
    a,b=contours[coarse],contours[fine];phia=arrays[coarse];phib=arrays[fine]
    restricted=phib.reshape(coarse,2,2*coarse,2).mean(axis=(1,3))
    stats[f'grid_{coarse}_{fine}']={
      'contour_vertex_hausdorff_m':max(directed(a,b),directed(b,a)),
      'phase_L1_over_domain':float(np.mean(np.abs(phia-restricted))),
      'note':'Contour vertex distances include sampling error; diffuse thickness changes with dx. No formal convergence order.'}
ax.set(xlabel='x [mm]',ylabel='y [mm]',title='Computed liquid fraction = 0.5 contours at 0.07 s\nGrid comparison; no paper shape reference',xlim=(2,8),ylim=(10,18))
ax.set_aspect('equal');ax.legend();fig.tight_layout();fig.savefig(OUT/'grid-shapes.png',dpi=180)
(OUT/'grid-field-review.json').write_text(json.dumps(stats,indent=2)+'\n')
print(json.dumps(stats,indent=2))
