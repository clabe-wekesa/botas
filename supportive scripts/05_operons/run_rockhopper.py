#!/usr/bin/env python3
from __future__ import annotations
import argparse,gzip,shutil,sys,tempfile
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'lib'))
from manifest import read_tsv,ensure_file
from common import require,run

def unpack(src,dst):
    if str(src).endswith('.gz'):
        with gzip.open(src,'rb') as fi, open(dst,'wb') as fo: shutil.copyfileobj(fi,fo)
    else: shutil.copyfile(src,dst)

def main():
    ap=argparse.ArgumentParser(description='Run Rockhopper on arbitrary paired-end samples. Manifest columns: sample,r1,r2,rockhopper_reference.')
    ap.add_argument('--manifest',required=True)
    ap.add_argument('--outdir',required=True)
    ap.add_argument('--jar',required=True)
    ap.add_argument('--threads',type=int,default=8)
    ap.add_argument('--java-mem',default='16g')
    ap.add_argument('--force',action='store_true')
    ap.add_argument('--dry-run',action='store_true')
    args=ap.parse_args(); require('java','/usr/bin/time'); jar=ensure_file(args.jar,'Rockhopper jar')
    rows=read_tsv(args.manifest,('sample','r1','r2','rockhopper_reference')); root=Path(args.outdir);root.mkdir(parents=True,exist_ok=True)
    for x in rows:
        s=x['sample'];r1=ensure_file(x['r1'],'R1');r2=ensure_file(x['r2'],'R2');refdir=Path(x['rockhopper_reference']).expanduser()
        if not refdir.is_dir(): raise FileNotFoundError(f'Rockhopper reference directory not found: {refdir}')
        for ext in ('fna','ptt','rnt'):
            if not (refdir/f'{s}.{ext}').is_file(): raise FileNotFoundError(refdir/f'{s}.{ext}')
        d=root/s;out=d/'results';out.mkdir(parents=True,exist_ok=True);operons=out/'_operons.txt'
        if operons.exists() and not args.force:
            print(f'[{s}] existing {operons}; skipping');continue
        if args.dry_run:
            print(f'[{s}] would unpack reads and run Rockhopper');continue
        with tempfile.TemporaryDirectory(prefix=f'rockhopper_{s}_') as td:
            tr1=Path(td)/f'{s}_1.fastq';tr2=Path(td)/f'{s}_2.fastq';unpack(r1,tr1);unpack(r2,tr2)
            cmd=['java',f'-Xmx{args.java_mem}','-cp',jar,'Rockhopper','-g',refdir,'-p',args.threads,'-y','true','-L',s,'-TIME','-o',out,f'{tr1}%{tr2}']
            with (d/'stdout.log').open('w') as so,(d/'stderr.log').open('w') as se:run(['/usr/bin/time','-v','-o',d/'time.txt',*cmd],stdout=so,stderr=se)
        if not operons.is_file(): raise RuntimeError(f'{s}: Rockhopper did not produce {operons}')
        print(f'[{s}] {operons}')
if __name__=='__main__':main()
