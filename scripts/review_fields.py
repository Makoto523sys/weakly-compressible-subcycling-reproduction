"""Plot saved fields at their actual snapshot times, without assigning target times."""
from pathlib import Path
import csv,json,tomllib
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

root=Path(__file__).resolve().parents[1]
out=root/'evidence/local-20260907'
def vtk(path,nx,ny):
    text=path.read_text()
    phi=np.fromstring(text.split('LOOKUP_TABLE default\n')[1].split('SCALARS')[0],sep=' ')
    p=np.fromstring(text.split('LOOKUP_TABLE default\n')[2].split('VECTORS')[0],sep=' ')
    v=np.fromstring(text.split('VECTORS velocity_m_per_s double\n')[1],sep=' ').reshape(-1,3)
    assert len(phi)==len(p)==len(v)==nx*ny
    return phi.reshape((nx,ny),order='F'),p.reshape((nx,ny),order='F'),v[:,:2].reshape((nx,ny,2),order='F')

fig,axes=plt.subplots(1,3,figsize=(9,7),sharex=True,sharey=True)
stats={}
for ax,name in zip(axes,['bubble-standard-mg','bubble-weak-initial','bubble-proposed-initial']):
    run=root/'results'/name;meta=tomllib.loads((run/'run.toml').read_text());c=meta['config']
    if (run/'snapshots.csv').exists():
        with (run/'snapshots.csv').open() as f: snapshots=list(csv.DictReader(f))
    else:
        with (run/'history.csv').open() as f: history=list(csv.DictReader(f))
        interval=meta['requested_t_end']/20;next_t=interval
        snapshots=[{'file':'field_000000.vtk','time':0.0}]
        for row in history[1:]:
            t=float(row['time'])
            if t>=next_t or t>=meta['requested_t_end']:
                snapshots.append({'file':f'field_{len(snapshots):06d}.vtk','time':t});next_t+=interval
        with (run/'snapshots-derived.csv').open('w') as f:
            writer=csv.DictWriter(f,fieldnames=['file','time']);writer.writeheader();writer.writerows(snapshots)
    last=snapshots[-1]
    # An actively written VTK can be incomplete: fall back to the preceding one.
    try:phi,p,v=vtk(run/last['file'],c['nx'],c['ny'])
    except (AssertionError,ValueError):
        last=snapshots[-2];phi,p,v=vtk(run/last['file'],c['nx'],c['ny'])
    x=(np.arange(c['nx'])+.5)*c['lx']/c['nx']*1000
    y=(np.arange(c['ny'])+.5)*c['ly']/c['ny']*1000
    # Fix the display scale; report raw overshoots below without changing data.
    ax.pcolormesh(x,y,phi.T,cmap='Blues',vmin=0,vmax=1,shading='nearest')
    ax.contour(x,y,phi.T,levels=[.5],colors='black',linewidths=.8)
    ax.set_aspect('equal');ax.set_title(f'{name}\nt={float(last["time"]):.6f} s',fontsize=9)
    ax.set_xlabel('x [mm]')
    stats[name]={'snapshot':last,'max_cell_speed_m_s':float(np.max(np.linalg.norm(v,axis=2))),
                 'pressure_min_Pa':float(p.min()),'pressure_max_Pa':float(p.max()),
                 'phi_min':float(phi.min()),'phi_max':float(phi.max())}
axes[0].set_ylabel('y [mm]')
fig.suptitle('Computed phase fields (color scale 0–1; raw extrema in JSON)\nActual times shown; no paper shape reference')
fig.tight_layout();fig.savefig(out/'bubble-fields.png',dpi=180)
(out/'field-review.json').write_text(json.dumps(stats,indent=2)+'\n')
