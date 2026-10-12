#!/usr/bin/env python3
import argparse,csv,subprocess,sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'lib'))
from manifest import read_tsv, ensure_file
from common import require, run

def main():
    ap=argparse.ArgumentParser(description='Run BOTAS getOperons on arbitrary BAM/GFF pairs.')
    ap.add_argument('--manifest',required=True,help='TSV columns: sample,bam,gff')
    ap.add_argument('--outdir',required=True)
    ap.add_argument('--threads',type=int,default=8)
    ap.add_argument('--gene-id',default='gene')
    ap.add_argument('--max-igd',type=int,default=70)
    ap.add_argument('--min-coverage',type=float,default=0.5)
    ap.add_argument('--min-cov-ratio',type=float,default=0.25)
    ap.add_argument('--min-score',type=float,default=0.0)
    ap.add_argument('--force',action='store_true'); ap.add_argument('--dry-run',action='store_true')
    args=ap.parse_args(); require('botas','/usr/bin/time')
    rows=read_tsv(args.manifest,('sample','bam','gff')); root=Path(args.outdir);root.mkdir(parents=True,exist_ok=True)
    outmanifest=root/'operon_outputs.tsv'
    with outmanifest.open('w',newline='') as fh:
        w=csv.writer(fh,delimiter='\t');w.writerow(['sample','operons'])
        for x in rows:
            s=x['sample'];bam=ensure_file(x['bam'],'BAM');gff=ensure_file(x['gff'],'GFF');d=root/s;d.mkdir(exist_ok=True);prefix=d/f'{s}_operons'
            expected=Path(str(prefix)+'.operons.botas')/'results'/f'{s}_operons.operons.tsv'
            if args.force or not expected.exists():
                cmd=['botas','getOperons','--bam',bam,'--gff',gff,'--threads',args.threads,'--gene-id',args.gene_id,'--strand-aware','--max-igd',args.max_igd,'--min-coverage',args.min_coverage,'--min-cov-ratio',args.min_cov_ratio,'--min-score',args.min_score,'--out',prefix]
                with (d/'stdout.log').open('w') as so,(d/'stderr.log').open('w') as se: run(['/usr/bin/time','-v','-o',d/'time.txt',*cmd],stdout=so,stderr=se,dry_run=args.dry_run)
            if not args.dry_run and not expected.exists(): raise RuntimeError(f'{s}: expected output missing: {expected}')
            w.writerow([s,expected])
    print(outmanifest)
if __name__=='__main__': main()
