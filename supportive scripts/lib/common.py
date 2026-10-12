from __future__ import annotations
import shlex, shutil, subprocess
from pathlib import Path

def require(*commands):
    missing=[c for c in commands if shutil.which(c) is None]
    if missing: raise RuntimeError('Missing required command(s): '+', '.join(missing))

def run(cmd, *, stdout=None, stderr=None, cwd=None, dry_run=False):
    print('+ ' + ' '.join(shlex.quote(str(x)) for x in cmd), flush=True)
    if dry_run: return
    subprocess.run([str(x) for x in cmd],check=True,stdout=stdout,stderr=stderr,cwd=cwd)

def fasta_ids(path):
    ids=[]
    opener=open
    import gzip
    if str(path).endswith('.gz'): opener=gzip.open
    with opener(path,'rt') as fh:
        for line in fh:
            if line.startswith('>'): ids.append(line[1:].split()[0])
    return ids

def ensure_parent(path): Path(path).parent.mkdir(parents=True,exist_ok=True)
