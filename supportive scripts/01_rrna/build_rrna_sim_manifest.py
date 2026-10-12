#!/usr/bin/env python3
import argparse,csv,sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'lib'))
from manifest import read_tsv

def main():
    ap=argparse.ArgumentParser(description='Build a sample manifest for rRNA simulated datasets created by build_rrna_benchmark.sh.')
    ap.add_argument('--source-manifest',required=True,help='Original sample manifest with sample,reference,gff')
    ap.add_argument('--benchmark-dir',required=True)
    ap.add_argument('--output',required=True)
    args=ap.parse_args(); root=Path(args.benchmark_dir); out=[]
    for r in read_tsv(args.source_manifest,('sample','reference','gff')):
        s=r['sample']; d=root/s
        r1=d/f'{s}_rrna_benchmark_100k_R1.fastq'; r2=d/f'{s}_rrna_benchmark_100k_R2.fastq'; truth=d/f'{s}_rrna_benchmark_truth.tsv'
        for p in (r1,r2,truth):
            if not p.is_file(): raise FileNotFoundError(p)
        out.append(dict(sample=s,r1=str(r1.resolve()),r2=str(r2.resolve()),reference=str(Path(r['reference']).resolve()),gff=str(Path(r['gff']).resolve()),truth=str(truth.resolve())))
    p=Path(args.output);p.parent.mkdir(parents=True,exist_ok=True)
    with p.open('w',newline='') as fh:w=csv.DictWriter(fh,fieldnames=['sample','r1','r2','reference','gff','truth'],delimiter='\t');w.writeheader();w.writerows(out)
    print(p)
if __name__=='__main__':main()
