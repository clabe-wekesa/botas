#!/usr/bin/env python3
import argparse,csv,sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'lib'))
from manifest import read_tsv

def main():
    ap=argparse.ArgumentParser(description='Create a BAM/GFF/reference manifest by selecting one aligner from an alignment summary and joining it to samples.tsv.')
    ap.add_argument('--samples',required=True,help='TSV columns: sample,gff,reference')
    ap.add_argument('--alignment-summary',required=True,help='TSV from run_4_aligners.py')
    ap.add_argument('--aligner',default='BOTAS')
    ap.add_argument('--output',required=True);args=ap.parse_args()
    meta={r['sample']:r for r in read_tsv(args.samples,('sample','gff','reference'))}; bams={}
    with open(args.alignment_summary,newline='') as fh:
        for r in csv.DictReader(fh,delimiter='\t'):
            if r.get('aligner')==args.aligner:bams[r['sample']]=r['bam']
    rows=[]
    for s,m in meta.items():
        if s not in bams: raise KeyError(f'No {args.aligner} BAM for {s}')
        rows.append(dict(sample=s,bam=str(Path(bams[s]).resolve()),gff=str(Path(m['gff']).resolve()),reference=str(Path(m['reference']).resolve())))
    p=Path(args.output);p.parent.mkdir(parents=True,exist_ok=True)
    with p.open('w',newline='') as fh:w=csv.DictWriter(fh,fieldnames=['sample','bam','gff','reference'],delimiter='\t');w.writeheader();w.writerows(rows)
    print(p)
if __name__=='__main__':main()
