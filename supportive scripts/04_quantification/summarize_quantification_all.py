#!/usr/bin/env python3
import argparse,csv,subprocess,sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'lib'))
from manifest import read_tsv,ensure_file

def count_fragments(bam):
    x=subprocess.check_output(['samtools','view','-c','-f','64','-F','2304',str(bam)],text=True).strip();return int(x)
def find_botas_tsv(d):
    c=[p for p in d.rglob('*.tsv') if 'summary' not in p.name.lower() and 'time' not in p.name.lower()]
    if not c: raise FileNotFoundError(f'No BOTAS result TSV under {d}')
    pref=[p for p in c if 'quant' in p.name.lower() or 'count' in p.name.lower()]
    return sorted(pref or c)[0]
def countcol(cols):
    lookup={c.lower():c for c in cols}
    for x in ['count','counts','fragment_count','fragments','assigned_fragments','raw_count','raw_counts']:
        if x in lookup:return lookup[x]
    raise ValueError('Cannot identify count column: '+','.join(cols))
def parse_botas(d,total):
    p=find_botas_tsv(d); assigned=0;features=0
    with p.open() as fh:
        r=csv.DictReader(fh,delimiter='\t'); col=countcol(r.fieldnames or [])
        for row in r:
            try:v=float(row[col])
            except:continue
            assigned+=v;features+=1
    assigned=round(assigned);return features,assigned,total-assigned
def parse_fc(p):
    features=sum(1 for x in open(p) if x.strip() and not x.startswith('#') and not x.startswith('Geneid'))
    s=Path(str(p)+'.summary');assigned=un=0
    with s.open() as fh:
        r=csv.reader(fh,delimiter='\t');next(r,None)
        for row in r:
            if len(row)<2:continue
            try:v=int(row[1])
            except:continue
            if row[0]=='Assigned':assigned+=v
            else:un+=v
    return features,assigned,un
def parse_ht(p):
    features=assigned=un=0
    with p.open() as fh:
        for line in fh:
            a=line.rstrip().split('\t')
            if len(a)!=2:continue
            try:v=int(a[1])
            except:continue
            if a[0].startswith('__'):un+=v
            else:features+=1;assigned+=v
    return features,assigned,un
def main():
    ap=argparse.ArgumentParser(description='Summarize outputs created by run_quantification.py.')
    ap.add_argument('--manifest',required=True,help='TSV columns: sample,bam,gff')
    ap.add_argument('--results-dir',required=True);ap.add_argument('--output',required=True);args=ap.parse_args()
    rows=[]; root=Path(args.results_dir)
    for x in read_tsv(args.manifest,('sample','bam','gff')):
        s=x['sample'];total=count_fragments(ensure_file(x['bam'],'BAM'));sd=root/s
        for tool in ['BOTAS','featureCounts','HTSeq']:
            try:
                if tool=='BOTAS':f,a,u=parse_botas(sd/'botas',total)
                elif tool=='featureCounts':f,a,u=parse_fc(sd/'featurecounts'/f'{s}.featureCounts.tsv')
                else:f,a,u=parse_ht(sd/'htseq'/f'{s}.htseq.tsv')
                den=a+u if tool!='BOTAS' else total
                rows.append(dict(sample=s,tool=tool,features=f,assigned=a,unassigned=u,total=den,assignment_pct=f'{100*a/den if den else 0:.2f}'))
            except FileNotFoundError: continue
    p=Path(args.output);p.parent.mkdir(parents=True,exist_ok=True)
    with p.open('w',newline='') as fh:w=csv.DictWriter(fh,fieldnames=['sample','tool','features','assigned','unassigned','total','assignment_pct'],delimiter='\t');w.writeheader();w.writerows(rows)
    print(p)
if __name__=='__main__':main()
