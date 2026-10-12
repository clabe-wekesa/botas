#!/usr/bin/env python3
import argparse,csv,sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'lib'))
from manifest import read_tsv,ensure_file

def attrs(s):
    d={}
    for x in s.strip().split(';'):
        if '=' in x:
            k,v=x.split('=',1);d[k]=v
    return d

def pair(a,b): return tuple(sorted((a,b)))
def protein_genes(gff):
    ids={};coding=set()
    with open(gff) as fh:
        for line in fh:
            if line.startswith('#'): continue
            f=line.rstrip().split('\t')
            if len(f)!=9: continue
            a=attrs(f[8])
            if f[2]=='gene':
                gid=a.get('ID');name=a.get('gene') or a.get('Name') or a.get('locus_tag') or gid
                if gid and name: ids[gid]=name
            elif f[2]=='CDS':
                p=a.get('Parent')
                if p: coding.update(p.split(','))
    return {ids[x] for x in coding if x in ids}
def load_pairs(path,column,valid):
    out=set()
    with open(path,newline='') as fh:
        for r in csv.DictReader(fh,delimiter='\t'):
            genes=[x.strip() for x in r[column].split(',') if x.strip()]
            for a,b in zip(genes,genes[1:]):
                if a in valid and b in valid: out.add(pair(a,b))
    return out
def main():
    ap=argparse.ArgumentParser(description='Compare BOTAS and Rockhopper adjacent-gene operon predictions for arbitrary samples.')
    ap.add_argument('--manifest',required=True,help='TSV columns: sample,gff,botas_operons,rockhopper_operons')
    ap.add_argument('--output',required=True);args=ap.parse_args()
    rows=[]
    for x in read_tsv(args.manifest,('sample','gff','botas_operons','rockhopper_operons')):
        valid=protein_genes(ensure_file(x['gff'],'GFF'))
        b=load_pairs(ensure_file(x['botas_operons'],'BOTAS operons'),'gene_ids',valid)
        r=load_pairs(ensure_file(x['rockhopper_operons'],'Rockhopper operons'),'Genes',valid)
        shared=b&r; union=b|r; p=len(shared)/len(b) if b else 0; rec=len(shared)/len(r) if r else 0; f1=2*p*rec/(p+rec) if p+rec else 0
        rows.append(dict(sample=x['sample'],botas_pairs=len(b),rockhopper_pairs=len(r),shared=len(shared),botas_only=len(b-r),rockhopper_only=len(r-b),botas_shared_pct=f'{100*p:.2f}',rockhopper_shared_pct=f'{100*rec:.2f}',agreement_f1=f'{f1:.4f}',jaccard=f'{len(shared)/len(union) if union else 0:.4f}'))
    p=Path(args.output);p.parent.mkdir(parents=True,exist_ok=True)
    with p.open('w',newline='') as fh:w=csv.DictWriter(fh,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
    print(p)
if __name__=='__main__':main()
