#!/usr/bin/env python3
from __future__ import annotations
import argparse, csv, subprocess, sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'lib'))
from manifest import read_tsv, ensure_file
from common import require, run

def timed(cmd,timefile,stdout=None,stderr=None,dry=False):
    run(['/usr/bin/time','-f','runtime_seconds\t%e\npeak_ram_kb\t%M','-o',timefile,*cmd],stdout=stdout,stderr=stderr,dry_run=dry)

def main():
    ap=argparse.ArgumentParser(description='Run BOTAS quantify, featureCounts and HTSeq on arbitrary BAM/GFF pairs.')
    ap.add_argument('--manifest',required=True,help='TSV columns: sample,bam,gff')
    ap.add_argument('--outdir',required=True)
    ap.add_argument('--threads',type=int,default=8)
    ap.add_argument('--feature-type',default='gene')
    ap.add_argument('--id-attribute',default='locus_tag')
    ap.add_argument('--min-overlap',type=int,default=1)
    ap.add_argument('--tools',nargs='+',choices=['botas','featurecounts','htseq'],default=['botas','featurecounts','htseq'])
    ap.add_argument('--force',action='store_true')
    ap.add_argument('--dry-run',action='store_true')
    args=ap.parse_args()
    needs=[]
    if 'botas' in args.tools: needs.append('botas')
    if 'featurecounts' in args.tools: needs.append('featureCounts')
    if 'htseq' in args.tools: needs.append('htseq-count')
    require(*needs,'/usr/bin/time')
    rows=read_tsv(args.manifest,('sample','bam','gff')); root=Path(args.outdir); root.mkdir(parents=True,exist_ok=True)
    manifest_out=root/'outputs.tsv'
    with manifest_out.open('w',newline='') as fh:
        w=csv.writer(fh,delimiter='\t');w.writerow(['sample','tool','output'])
        for x in rows:
            s=x['sample']; bam=ensure_file(x['bam'],'BAM'); gff=ensure_file(x['gff'],'GFF'); sd=root/s; sd.mkdir(exist_ok=True)
            if 'botas' in args.tools:
                d=sd/'botas';d.mkdir(exist_ok=True);prefix=d/f'{s}_quant'; marker=d/'completed.txt'
                if args.force or not marker.exists():
                    with (d/'stdout.log').open('w') as so,(d/'stderr.log').open('w') as se:
                        timed(['botas','quantify','-b',bam,'-g',gff,'--feature-type',args.feature_type,'--id-attribute',args.id_attribute,'--min-overlap',args.min_overlap,'-o',prefix],d/'time.tsv',so,se,args.dry_run)
                    if not args.dry_run: marker.write_text('ok\n')
                w.writerow([s,'BOTAS',prefix])
            if 'featurecounts' in args.tools:
                d=sd/'featurecounts';d.mkdir(exist_ok=True);out=d/f'{s}.featureCounts.tsv'
                if args.force or not out.exists():
                    with (d/'stdout.log').open('w') as so,(d/'stderr.log').open('w') as se:
                        timed(['featureCounts','-T',args.threads,'-p','--countReadPairs','-t',args.feature_type,'-g',args.id_attribute,'--minOverlap',args.min_overlap,'-s','0','-a',gff,'-o',out,bam],d/'time.tsv',so,se,args.dry_run)
                w.writerow([s,'featureCounts',out])
            if 'htseq' in args.tools:
                d=sd/'htseq';d.mkdir(exist_ok=True);out=d/f'{s}.htseq.tsv'
                if args.force or not out.exists():
                    with out.open('w') as so,(d/'stderr.log').open('w') as se:
                        timed(['htseq-count','-f','bam','-r','pos','-s','no','-t',args.feature_type,'-i',args.id_attribute,'-m','union',bam,gff],d/'time.tsv',so,se,args.dry_run)
                w.writerow([s,'HTSeq',out])
    print(manifest_out)
if __name__=='__main__': main()
