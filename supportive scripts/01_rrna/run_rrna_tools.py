#!/usr/bin/env python3
from __future__ import annotations
import argparse,csv,sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'lib'))
from manifest import read_tsv, ensure_file
from common import require, run, fasta_ids

def timed(cmd,timefile,stdout=None,stderr=None,dry=False):
    run(['/usr/bin/time','-f','runtime_seconds\t%e\npeak_ram_kb\t%M','-o',timefile,*cmd],stdout=stdout,stderr=stderr,dry_run=dry)

def botas_bam(prefix):
    r=Path(str(prefix)+'.botas')/'results'
    xs=sorted(r.glob('*.bam'))
    return xs[0] if xs else None

def first_existing(candidates):
    return next((p for p in candidates if p.exists()),None)

def main():
    ap=argparse.ArgumentParser(description='Run an rRNA filtering benchmark tool on arbitrary paired-end samples from a TSV manifest.')
    ap.add_argument('--tool',required=True,choices=['botas','sortmerna','ribodetector'])
    ap.add_argument('--manifest',required=True,help='TSV: sample,r1,r2; BOTAS also requires reference')
    ap.add_argument('--outdir',required=True)
    ap.add_argument('--threads',type=int,default=8)
    ap.add_argument('--rrna-ref',help='SortMeRNA reference FASTA')
    ap.add_argument('--rrna-k',type=int,default=18)
    ap.add_argument('--rrna-min-hits',type=int,default=10)
    ap.add_argument('--read-length',type=int,default=150)
    ap.add_argument('--circular',choices=['all','none'],default='all')
    ap.add_argument('--force',action='store_true')
    ap.add_argument('--dry-run',action='store_true')
    args=ap.parse_args()
    if args.tool=='botas': require('botas','/usr/bin/time'); required=('sample','r1','r2','reference')
    elif args.tool=='sortmerna': require('sortmerna','/usr/bin/time'); required=('sample','r1','r2')
    else: require('ribodetector_cpu','/usr/bin/time'); required=('sample','r1','r2')
    if args.tool=='sortmerna' and not args.rrna_ref: ap.error('--rrna-ref is required for SortMeRNA')
    rows=read_tsv(args.manifest,required); root=Path(args.outdir); root.mkdir(parents=True,exist_ok=True); results=[]
    for x in rows:
        s=x['sample']; r1=ensure_file(x['r1'],'R1'); r2=ensure_file(x['r2'],'R2'); d=root/s/args.tool; d.mkdir(parents=True,exist_ok=True)
        retained_bam=rrna_r1=non_r1=''
        if args.tool=='botas':
            ref=ensure_file(x['reference'],'reference'); contigs=','.join(fasta_ids(ref)) if args.circular=='all' else ''
            idxp=d/f'{s}.idx'; idx=Path(str(idxp)+'.botas')/'results'/f'{s}.idx.botas.idx'
            if args.force or not idx.exists():
                cmd=['botas','index','-r',ref,'-o',idxp]
                if contigs: cmd += ['--circular-contigs',contigs]
                timed(cmd,d/'index.time.tsv',dry=args.dry_run)
            outp=d/f'{s}.filtered'; bam=botas_bam(outp)
            if args.force or bam is None:
                cmd=['botas','align','--ref',ref,'--index',idx,'--fq1',r1,'--fq2',r2,'-t',args.threads,'--pool','--filter-rrna','--rrna-k',args.rrna_k,'--rrna-min-hits',args.rrna_min_hits,'-o',outp]
                if contigs: cmd += ['--circular-contigs',contigs]
                with (d/'stdout.log').open('w') as so,(d/'stderr.log').open('w') as se: timed(cmd,d/'time.tsv',so,se,args.dry_run)
                bam=botas_bam(outp)
            if not args.dry_run and bam is None: raise RuntimeError(f'{s}: BOTAS retained-read BAM not found')
            retained_bam=str(bam.resolve()) if bam else ''
        elif args.tool=='sortmerna':
            rrna=ensure_file(args.rrna_ref,'rRNA reference'); work=d/'work'; work.mkdir(exist_ok=True)
            cmd=['sortmerna','--ref',rrna,'--reads',r1,'--reads',r2,'--threads',args.threads,'--workdir',work,'--fastx','--paired_in','--out2','--aligned',d/'rrna','--other',d/'non_rrna']
            with (d/'stdout.log').open('w') as so,(d/'stderr.log').open('w') as se: timed(cmd,d/'time.tsv',so,se,args.dry_run)
            if not args.dry_run:
                rrnafile=first_existing([d/'rrna_fwd.fq',d/'rrna_fwd.fastq',d/'rrna_1.fq',d/'rrna_1.fastq'])
                nonfile=first_existing([d/'non_rrna_fwd.fq',d/'non_rrna_fwd.fastq',d/'non_rrna_1.fq',d/'non_rrna_1.fastq'])
                if rrnafile is None or nonfile is None: raise RuntimeError(f'{s}: could not locate SortMeRNA paired output FASTQ files in {d}')
                rrna_r1=str(rrnafile.resolve()); non_r1=str(nonfile.resolve())
        else:
            rr1=d/'rrna_1.fastq'; rr2=d/'rrna_2.fastq'; nr1=d/'nonrrna_1.fastq'; nr2=d/'nonrrna_2.fastq'
            cmd=['ribodetector_cpu','-l',args.read_length,'-i',r1,r2,'-o',nr1,nr2,'-r',rr1,rr2,'-e','none','-t',args.threads,'--log',d/'ribodetector.log']
            timed(cmd,d/'time.tsv',dry=args.dry_run)
            if not args.dry_run:
                rrna_r1=str(rr1.resolve());non_r1=str(nr1.resolve())
        results.append(dict(sample=s,tool=args.tool,retained_bam=retained_bam,rrna_r1=rrna_r1,non_rrna_r1=non_r1,output_dir=str(d.resolve())))
    mf=root/f'{args.tool}_outputs.tsv'
    with mf.open('w',newline='') as fh:
        w=csv.DictWriter(fh,fieldnames=['sample','tool','retained_bam','rrna_r1','non_rrna_r1','output_dir'],delimiter='\t');w.writeheader();w.writerows(results)
    print(mf)
if __name__=='__main__': main()
