#!/usr/bin/env python3
import argparse,csv,sys
from collections import defaultdict
from pathlib import Path
import pysam
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'lib'))
from manifest import read_tsv,ensure_file

def usable(r): return (not r.is_unmapped and r.is_paired and not r.is_secondary and not r.is_supplementary and not r.is_qcfail)
def main():
    ap=argparse.ArgumentParser(description='Count unique paired-end fragments with at least one primary XC:Z:JUNCTION alignment.')
    ap.add_argument('--manifest',required=True,help='TSV columns: sample,bam')
    ap.add_argument('--output',required=True);ap.add_argument('--details');ap.add_argument('--min-mapq',type=int,default=20);args=ap.parse_args()
    rows=[];detail=[]
    for x in read_tsv(args.manifest,('sample','bam')):
        sample=x['sample']; bam=ensure_file(x['bam'],'BAM'); frag=defaultdict(list); examined=set()
        with pysam.AlignmentFile(str(bam),'rb') as bf:
            for r in bf.fetch(until_eof=True):
                if not usable(r): continue
                examined.add(r.query_name)
                if r.has_tag('XC') and r.get_tag('XC')=='JUNCTION':
                    frag[r.query_name].append(r)
        high=sum(1 for rs in frag.values() if max(r.mapping_quality for r in rs)>=args.min_mapq)
        one=sum(1 for rs in frag.values() if len(rs)==1); two=sum(1 for rs in frag.values() if len(rs)>=2)
        rows.append(dict(sample=sample,bam=str(bam),examined_fragments=len(examined),junction_fragments=len(frag),junction_fragments_mapq_ge_min=high,one_tagged_mate=one,two_or_more_tagged_mates=two,min_mapq=args.min_mapq))
        if args.details:
            for q,rs in sorted(frag.items()):
                detail.append(dict(sample=sample,fragment=q,tagged_alignments=len(rs),max_mapq=max(r.mapping_quality for r in rs),references=','.join(sorted({r.reference_name for r in rs if r.reference_name}))))
    p=Path(args.output);p.parent.mkdir(parents=True,exist_ok=True)
    with p.open('w',newline='') as fh:w=csv.DictWriter(fh,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
    if args.details:
        d=Path(args.details);d.parent.mkdir(parents=True,exist_ok=True)
        with d.open('w',newline='') as fh:w=csv.DictWriter(fh,fieldnames=['sample','fragment','tagged_alignments','max_mapq','references'],delimiter='\t');w.writeheader();w.writerows(detail)
    print(p)
if __name__=='__main__':main()
