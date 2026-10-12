#!/usr/bin/env python3
import argparse,csv,gzip,subprocess,sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'lib'))
from manifest import read_tsv,ensure_file

def fastq_records(path):
    p=Path(path); op=gzip.open if str(p).endswith('.gz') else open; n=0
    with op(p,'rt') as fh:
        for _ in fh: n+=1
    if n%4: raise RuntimeError(f'FASTQ line count not divisible by 4: {p}')
    return n//4

def bam_pairs(path):
    out=subprocess.check_output(['samtools','view','-c','-f','64','-F','2304',str(path)],text=True)
    return int(out.strip())

def main():
    ap=argparse.ArgumentParser(description='Summarize experimental rRNA filtering from run_rrna_tools.py output manifests.')
    ap.add_argument('--samples',required=True,help='Original sample TSV columns: sample,r1,r2')
    ap.add_argument('--predictions-manifest',required=True,help='Output TSV from run_rrna_tools.py')
    ap.add_argument('--output',required=True);args=ap.parse_args()
    samples={x['sample']:x for x in read_tsv(args.samples,('sample','r1','r2'))};rows=[]
    for pred in read_tsv(args.predictions_manifest,('sample','tool')):
        s=pred['sample']; tool=pred['tool']
        if s not in samples: raise KeyError(f'No input sample metadata for {s}')
        input_pairs=fastq_records(ensure_file(samples[s]['r1'],'R1'))
        if pred.get('retained_bam'):
            retained=bam_pairs(ensure_file(pred['retained_bam'],'retained BAM')); filtered=input_pairs-retained
        else:
            filtered=fastq_records(ensure_file(pred['rrna_r1'],'rRNA FASTQ')); retained=fastq_records(ensure_file(pred['non_rrna_r1'],'non-rRNA FASTQ'))
            if filtered+retained!=input_pairs:
                raise RuntimeError(f'{s}: classified pairs ({filtered+retained}) != input pairs ({input_pairs})')
        rows.append(dict(sample=s,tool=tool,input_pairs=input_pairs,filtered_pairs=filtered,filtered_pct=f'{100*filtered/input_pairs if input_pairs else 0:.3f}',retained_pairs=retained,retained_pct=f'{100*retained/input_pairs if input_pairs else 0:.3f}'))
    p=Path(args.output);p.parent.mkdir(parents=True,exist_ok=True)
    with p.open('w',newline='') as fh:
        w=csv.DictWriter(fh,fieldnames=['sample','tool','input_pairs','filtered_pairs','filtered_pct','retained_pairs','retained_pct'],delimiter='\t');w.writeheader();w.writerows(rows)
    print(p)
if __name__=='__main__':main()
