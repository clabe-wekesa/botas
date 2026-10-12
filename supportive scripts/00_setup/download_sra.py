#!/usr/bin/env python3
import argparse, shutil, subprocess, tempfile
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'lib'))
from manifest import read_tsv

def run(cmd):
    subprocess.run(cmd, check=True)

def main():
    ap=argparse.ArgumentParser(description='Download paired-end SRA runs from a two-column sample/sra TSV manifest.')
    ap.add_argument('--manifest', required=True, help='TSV with columns: sample, sra')
    ap.add_argument('--outdir', required=True)
    ap.add_argument('--threads', type=int, default=8)
    ap.add_argument('--keep-fastq', action='store_true', help='Do not gzip FASTQ outputs')
    args=ap.parse_args()
    if shutil.which('fasterq-dump') is None: ap.error('fasterq-dump not found')
    compressor='pigz' if shutil.which('pigz') else 'gzip'
    rows=read_tsv(args.manifest, ('sample','sra'))
    out=Path(args.outdir); out.mkdir(parents=True,exist_ok=True)
    for row in rows:
        sample,runid=row['sample'],row['sra']
        print(f'[{sample}] {runid}', flush=True)
        with tempfile.TemporaryDirectory(prefix=f'{runid}.',dir=out) as td:
            run(['fasterq-dump',runid,'--split-files','--threads',str(args.threads),'--outdir',td])
            r1=Path(td)/f'{runid}_1.fastq'; r2=Path(td)/f'{runid}_2.fastq'
            if not (r1.is_file() and r2.is_file()): raise RuntimeError(f'{runid}: paired FASTQ files were not produced')
            dst1=out/f'{sample}_1.fastq'; dst2=out/f'{sample}_2.fastq'
            shutil.move(r1,dst1); shutil.move(r2,dst2)
            if not args.keep_fastq:
                if compressor=='pigz':
                    run(['pigz','-f','-p',str(args.threads),str(dst1),str(dst2)])
                else:
                    run(['gzip','-f',str(dst1),str(dst2)])
    print(f'Finished {len(rows)} run(s). Output: {out}')
if __name__=='__main__': main()
