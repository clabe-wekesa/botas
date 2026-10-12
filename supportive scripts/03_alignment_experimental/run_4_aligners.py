#!/usr/bin/env python3
from __future__ import annotations
import argparse, csv, os, shlex, subprocess, sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'lib'))
from manifest import read_tsv, ensure_file
from common import require, run, fasta_ids

def timed(cmd,timefile,stdout=None,stderr=None,dry=False):
    full=['/usr/bin/time','-f','runtime_seconds\t%e\npeak_ram_kb\t%M','-o',str(timefile),*map(str,cmd)]
    run(full,stdout=stdout,stderr=stderr,dry_run=dry)

def find_botas_bam(prefix):
    p=Path(str(prefix)+'.botas')/'results'
    bams=list(p.glob('*.bam'))
    return bams[0] if bams else None

def main():
    ap=argparse.ArgumentParser(description='Run BOTAS, Bowtie2, BWA-MEM2 and minimap2 on arbitrary paired-end samples from a TSV manifest.')
    ap.add_argument('--manifest',required=True,help='TSV columns: sample,r1,r2,reference (gff may also be present)')
    ap.add_argument('--outdir',required=True)
    ap.add_argument('--threads',type=int,default=8)
    ap.add_argument('--circular',choices=['all','none'],default='all',help='Treat all FASTA records as circular for BOTAS, or none.')
    ap.add_argument('--botas-overhang-percent',type=int,default=5)
    ap.add_argument('--force',action='store_true')
    ap.add_argument('--dry-run',action='store_true')
    args=ap.parse_args()
    require('botas','bowtie2','bowtie2-build','bwa-mem2','minimap2','samtools','/usr/bin/time')
    rows=read_tsv(args.manifest,('sample','r1','r2','reference'))
    root=Path(args.outdir); root.mkdir(parents=True,exist_ok=True)
    summary=root/'alignment_benchmark_summary.tsv'
    with summary.open('w',newline='') as sf:
        w=csv.writer(sf,delimiter='\t'); w.writerow(['sample','aligner','bam','flagstat'])
        for row in rows:
            sample=row['sample']; r1=ensure_file(row['r1'],'R1'); r2=ensure_file(row['r2'],'R2'); ref=ensure_file(row['reference'],'reference')
            sdir=root/sample; sdir.mkdir(parents=True,exist_ok=True)
            contigs=','.join(fasta_ids(ref)) if args.circular=='all' else ''
            # BOTAS
            bdir=sdir/'botas'; bdir.mkdir(exist_ok=True)
            idxprefix=bdir/f'{sample}.idx'; idx=Path(str(idxprefix)+'.botas')/'results'/f'{sample}.idx.botas.idx'
            if args.force or not idx.exists():
                cmd=['botas','index','-r',ref,'-o',idxprefix]
                if contigs: cmd += ['--circular-contigs',contigs,'--circular-overhang-percent',args.botas_overhang_percent]
                timed(cmd,bdir/'index.time.tsv',dry=args.dry_run)
            outprefix=bdir/sample
            bam=find_botas_bam(outprefix)
            if args.force or bam is None:
                cmd=['botas','align','--ref',ref,'--index',idx,'--fq1',r1,'--fq2',r2,'-t',args.threads,'--pool','-o',outprefix]
                if contigs: cmd += ['--circular-contigs',contigs]
                timed(cmd,bdir/'align.time.tsv',dry=args.dry_run)
                bam=find_botas_bam(outprefix)
            if not args.dry_run and bam is None: raise RuntimeError(f'{sample}: BOTAS BAM not found')
            if bam:
                sorted_bam=bdir/f'{sample}.sorted.bam'
                if args.force or not sorted_bam.exists(): run(['samtools','sort','-@',args.threads,'-o',sorted_bam,bam],dry_run=args.dry_run)
                if not args.dry_run: run(['samtools','index','-@',args.threads,sorted_bam]); fs=bdir/'flagstat.txt'; fs.write_text(subprocess.check_output(['samtools','flagstat',str(sorted_bam)],text=True)); w.writerow([sample,'BOTAS',sorted_bam,fs])
            # Bowtie2
            d=sdir/'bowtie2'; d.mkdir(exist_ok=True); pref=d/'index'; bam=d/f'{sample}.bam'
            if args.force or not (d/'index.1.bt2').exists() and not (d/'index.1.bt2l').exists(): run(['bowtie2-build',ref,pref],dry_run=args.dry_run)
            if args.force or not bam.exists():
                shell=f"bowtie2 -x {shlex.quote(str(pref))} -1 {shlex.quote(str(r1))} -2 {shlex.quote(str(r2))} -p {args.threads} | samtools view -@ {args.threads} -b -o {shlex.quote(str(bam))}"
                timed(['bash','-o','pipefail','-c',shell],d/'align.time.tsv',dry=args.dry_run)
            if not args.dry_run and bam.exists(): fs=d/'flagstat.txt'; fs.write_text(subprocess.check_output(['samtools','flagstat',str(bam)],text=True)); w.writerow([sample,'Bowtie2',bam,fs])
            # BWA-MEM2
            d=sdir/'bwa_mem2'; d.mkdir(exist_ok=True); bref=d/'reference.fna'; bam=d/f'{sample}.bam'
            if not bref.exists() and not args.dry_run: bref.write_bytes(ref.read_bytes())
            if args.force or not any(d.glob('reference.fna.*')): run(['bwa-mem2','index',bref],dry_run=args.dry_run)
            if args.force or not bam.exists():
                shell=f"bwa-mem2 mem -t {args.threads} {shlex.quote(str(bref))} {shlex.quote(str(r1))} {shlex.quote(str(r2))} | samtools view -@ {args.threads} -b -o {shlex.quote(str(bam))}"
                timed(['bash','-o','pipefail','-c',shell],d/'align.time.tsv',dry=args.dry_run)
            if not args.dry_run and bam.exists(): fs=d/'flagstat.txt'; fs.write_text(subprocess.check_output(['samtools','flagstat',str(bam)],text=True)); w.writerow([sample,'BWA-MEM2',bam,fs])
            # minimap2
            d=sdir/'minimap2'; d.mkdir(exist_ok=True); mmi=d/'reference.mmi'; bam=d/f'{sample}.bam'
            if args.force or not mmi.exists(): run(['minimap2','-d',mmi,ref],dry_run=args.dry_run)
            if args.force or not bam.exists():
                shell=f"minimap2 -ax sr -t {args.threads} {shlex.quote(str(mmi))} {shlex.quote(str(r1))} {shlex.quote(str(r2))} | samtools view -@ {args.threads} -b -o {shlex.quote(str(bam))}"
                timed(['bash','-o','pipefail','-c',shell],d/'align.time.tsv',dry=args.dry_run)
            if not args.dry_run and bam.exists(): fs=d/'flagstat.txt'; fs.write_text(subprocess.check_output(['samtools','flagstat',str(bam)],text=True)); w.writerow([sample,'Minimap2',bam,fs])
    print(summary)
if __name__=='__main__': main()
