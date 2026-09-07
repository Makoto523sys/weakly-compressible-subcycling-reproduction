"""Finish only the local evidence collection after the already-started refinement.

This does not approve reproduction, start porous work, or publish anything.
"""
from pathlib import Path
import os,subprocess,sys,time,tomllib
root=Path(__file__).resolve().parents[1]
meta=root/'results/bubble-proposed-512-fixed/run.toml'
while True:
    try:status=tomllib.loads(meta.read_text())['status']
    except (FileNotFoundError,tomllib.TOMLDecodeError):
        time.sleep(30);continue
    if status!='running_unvalidated':break
    time.sleep(30)
print('Refinement status:',status,flush=True)
if status!='completed_unvalidated':sys.exit(1)
env=dict(os.environ,MPLCONFIGDIR=str(root/'.mplconfig'))
for script in ('review_sensitivity.py','review_local.py','review_fields.py','review_grid_fields.py','collect_evidence.py'):
    subprocess.run([sys.executable,str(root/'scripts'/script)],cwd=root,env=env,check=True)
print('Refinement evidence updated; scientific review and porous implementation remain pending.',flush=True)
