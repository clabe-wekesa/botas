#!/usr/bin/env python3
import argparse,csv,sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'lib'))
from manifest import read_tsv, ensure_file

def attrs(s):
    d={}
    for x in s.strip().split(';'):
        if '=' in x:
            k,v=x.split('=',1);d[k]=v
    return d

def pair(a,b): return tuple(sorted((a,b)))
def protein_names(gff):
    ids={};coding=set()
    with open(gff) as fh:
        for line in fh:
            if line.startswith('#'): continue
            f=line.rstrip().split('\t')
            if len(f)!=9: continue
            a=attrs(f[8])
            if f[2]=='gene':
                gid=a.get('ID'); name=a.get('gene') or a.get('Name') or a.get('locus_tag') or gid
                if gid and name: ids[gid]=name
            elif f[2]=='CDS':
                par=a.get('Parent')
                if par: coding.update(par.split(','))
    return {ids[x] for x in coding if x in ids}
def truth_pairs(path,valid):
    out=set()
    with open(path,newline='') as fh:
        for r in csv.DictReader(fh,delimiter='\t'):
            a,b=r['gene1'].strip(),r['gene2'].strip()
            if a in valid and b in valid: out.add(pair(a,b))
    return out
def botas_pairs(path,valid):
    out=set()
    with open(path,newline='') as fh:
        for r in csv.DictReader(fh,delimiter='\t'):
            genes=[x.strip() for x in r['gene_ids'].split(',') if x.strip()]
            for a,b in zip(genes,genes[1:]):
                if a in valid and b in valid: out.add(pair(a,b))
    return out
def rockhopper_pairs(path,valid):
    out=set()
    with open(path,newline='') as fh:
        for r in csv.DictReader(fh,delimiter='\t'):
            genes=[x.strip() for x in r['Genes'].split(',') if x.strip()]
            for a,b in zip(genes,genes[1:]):
                if a in valid and b in valid: out.add(pair(a,b))
    return out
def main():
    ap=argparse.ArgumentParser(description='Build curated-operon benchmark summary for any samples.')
    ap.add_argument('--manifest',required=True,help='TSV columns: sample,reference_name,gff,truth_pairs,botas_operons,rockhopper_operons')
    ap.add_argument('--output',required=True);args=ap.parse_args()
    rows=read_tsv(args.manifest,('sample','reference_name','gff','truth_pairs','botas_operons','rockhopper_operons'));out=[]
    for x in rows:
        valid=protein_names(ensure_file(x['gff'],'GFF'));truth=truth_pairs(ensure_file(x['truth_pairs'],'truth pairs'),valid)
        for method,pred in [('BOTAS',botas_pairs(ensure_file(x['botas_operons'],'BOTAS operons'),valid)),('Rockhopper',rockhopper_pairs(ensure_file(x['rockhopper_operons'],'Rockhopper operons'),valid))]:
            sup=pred & truth; absent=pred-truth; missed=truth-pred
            out.append(dict(sample=x['sample'],reference=x['reference_name'],method=method,evaluable_curated_pairs=len(truth),predicted_pairs=len(pred),curated_supported_pairs=len(sup),predictions_not_represented=len(absent),missed_curated_pairs=len(missed),support_fraction_pct=f'{100*len(sup)/len(pred) if pred else 0:.2f}',curated_recall_pct=f'{100*len(sup)/len(truth) if truth else 0:.2f}'))
    p=Path(args.output);p.parent.mkdir(parents=True,exist_ok=True)
    with p.open('w',newline='') as fh:w=csv.DictWriter(fh,fieldnames=list(out[0]),delimiter='\t');w.writeheader();w.writerows(out)
    print(p)
if __name__=='__main__':main()
