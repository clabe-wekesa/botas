#!/usr/bin/env python3
import argparse,csv
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'lib'))
from manifest import read_tsv

def pick(base, patterns):
    for pat in patterns:
        p=base/pat
        if p.exists(): return str(p.resolve())
    return ''

def main():
    ap=argparse.ArgumentParser(description='Build a sample manifest from sample names plus standard FASTQ/reference locations.')
    ap.add_argument('--samples', required=True, help='TSV containing at least a sample column')
    ap.add_argument('--fastq-dir', required=True)
    ap.add_argument('--reference-dir', required=True)
    ap.add_argument('--output', required=True)
    args=ap.parse_args()
    rows=read_tsv(args.samples,('sample',)); fq=Path(args.fastq_dir); rr=Path(args.reference_dir)
    out=[]
    for x in rows:
        s=x['sample']; d=rr/s
        r1=pick(fq,[f'{s}_1.fastq.gz',f'{s}_R1.fastq.gz',f'{s}_1.fastq',f'{s}_R1.fastq'])
        r2=pick(fq,[f'{s}_2.fastq.gz',f'{s}_R2.fastq.gz',f'{s}_2.fastq',f'{s}_R2.fastq'])
        ref=pick(d,[f'{s}.fna',f'{s}.fa',f'{s}.fasta'])
        gff=pick(d,[f'{s}.gff3',f'{s}.gff'])
        if not all([r1,r2,ref,gff]): raise SystemExit(f'{s}: could not resolve r1/r2/reference/gff')
        out.append(dict(sample=s,r1=r1,r2=r2,reference=ref,gff=gff))
    p=Path(args.output); p.parent.mkdir(parents=True,exist_ok=True)
    with p.open('w',newline='') as fh:
        w=csv.DictWriter(fh,fieldnames=['sample','r1','r2','reference','gff'],delimiter='\t');w.writeheader();w.writerows(out)
    print(p)
if __name__=='__main__': main()
