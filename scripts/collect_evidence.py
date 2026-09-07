"""Archive finished histories and verify their solver hashes; omit large VTK files."""
from pathlib import Path
import hashlib, shutil, tomllib

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'evidence/local-20260907'
def digest(src):
    h=hashlib.sha256()
    for p in sorted(src.rglob('*.jl')):
        h.update(str(p.relative_to(src.parent)).encode()+b'\0'+p.read_bytes()+b'\0')
    return h.hexdigest()

known={}
for p in sorted((ROOT/'results').glob('*/src')):
    known[digest(p)]=p

for run in sorted((ROOT/'results').iterdir()):
    if not run.is_dir() or (run/'PROVENANCE_WARNING.md').exists():continue
    info=run/'run.toml'
    if not info.exists():info=run/'summary.toml'
    if not info.exists():continue
    meta=tomllib.loads(info.read_text())
    if meta['status'] not in ('completed_unvalidated','measured_requires_review',
                               'interrupted_after_multigrid_equivalence_check'):continue
    hashes=[meta['source_sha256']] if 'source_sha256' in meta else [r['source_sha256'] for r in meta['runs']]
    assert len(set(hashes))==1, f'Mixed revisions: {run}'
    src=known[hashes[0]]
    target=OUT/'runs'/run.name;target.mkdir(parents=True,exist_ok=True)
    for p in run.iterdir():
        if p.suffix in ('.csv','.toml'):shutil.copy2(p,target/p.name)
    shutil.copytree(src,target/'src',dirs_exist_ok=True)
    if not (target/'Project.toml').exists() and (src.parent/'Project.toml').exists():
        shutil.copy2(src.parent/'Project.toml',target/'Project.toml')
    assert digest(target/'src')==hashes[0]
    print(f'{run.name}: {meta["status"]}, source verified {hashes[0]}')
