#!/usr/bin/env python3
import argparse,csv,gzip,math,subprocess,sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'lib'))
from manifest import read_tsv,ensure_file

def norm_id(x):
    x=x.split()[0].lstrip('@')
    if x.endswith('/1') or x.endswith('/2'): x=x[:-2]
    return x

def read_truth(path):
    with open(path,newline='') as fh:
        r=csv.DictReader(fh,delimiter='\t')
        req={'pair_id','truth_class'}
        if not req.issubset(r.fieldnames or []): raise ValueError(f'{path}: truth TSV needs pair_id and truth_class')
        return {row['pair_id']:row['truth_class'] for row in r}

def fastq_ids(path):
    op=gzip.open if str(path).endswith('.gz') else open; ids=set()
    with op(path,'rt') as f:
        while True:
            h=f.readline()
            if not h: break
            seq=f.readline();plus=f.readline();qual=f.readline()
            if not qual: raise RuntimeError(f'incomplete FASTQ: {path}')
            ids.add(norm_id(h))
    return ids

def bam_ids(path):
    p=subprocess.Popen(['samtools','view',str(path)],stdout=subprocess.PIPE,text=True)
    ids={norm_id(line.split('\t',1)[0]) for line in p.stdout if line.strip()}
    if p.wait()!=0: raise RuntimeError('samtools view failed')
    return ids

def calc(tp,fp,tn,fn):
    sens=tp/(tp+fn) if tp+fn else 0; spec=tn/(tn+fp) if tn+fp else 0; prec=tp/(tp+fp) if tp+fp else 0; rec=sens
    f1=2*prec*rec/(prec+rec) if prec+rec else 0; acc=(tp+tn)/(tp+tn+fp+fn) if tp+tn+fp+fn else 0
    den=math.sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn));mcc=(tp*tn-fp*fn)/den if den else 0
    return dict(TP=tp,FP=fp,TN=tn,FN=fn,sensitivity=sens,specificity=spec,precision=prec,recall=rec,F1=f1,accuracy=acc,false_positive_rate=(fp/(fp+tn) if fp+tn else 0),MCC=mcc)

def score(truth_path,pred):
    truth=read_truth(truth_path); tids=set(truth)
    if pred.get('retained_bam'):
        non=bam_ids(ensure_file(pred['retained_bam'],'retained BAM')); rr=tids-non
    else:
        rr=fastq_ids(ensure_file(pred['rrna_r1'],'rRNA FASTQ')); non=fastq_ids(ensure_file(pred['non_rrna_r1'],'non-rRNA FASTQ'))
        if rr & non: raise RuntimeError('IDs occur in both predicted classes')
        missing=tids-(rr|non);extra=(rr|non)-tids
        if missing or extra: raise RuntimeError(f'classification/truth ID mismatch: missing={len(missing)} extra={len(extra)}')
    extra=(rr|non)-tids
    if extra: raise RuntimeError(f'{len(extra)} predicted IDs absent from truth')
    return calc(sum(truth[x]=='rRNA' for x in rr),sum(truth[x]=='non_rRNA' for x in rr),sum(truth[x]=='non_rRNA' for x in non),sum(truth[x]=='rRNA' for x in non))

def main():
    ap=argparse.ArgumentParser(description='Score rRNA predictions against simulated truth, either for one sample or in batch.')
    ap.add_argument('--truth');ap.add_argument('--retained-bam');ap.add_argument('--rrna-fastq');ap.add_argument('--non-rrna-fastq');ap.add_argument('--sample',default='sample')
    ap.add_argument('--truth-manifest',help='TSV columns: sample,truth')
    ap.add_argument('--predictions-manifest',help='Output TSV from run_rrna_tools.py')
    ap.add_argument('--output',required=True);args=ap.parse_args(); rows=[]
    if args.truth_manifest or args.predictions_manifest:
        if not (args.truth_manifest and args.predictions_manifest): ap.error('batch mode requires --truth-manifest and --predictions-manifest')
        truth={x['sample']:x['truth'] for x in read_tsv(args.truth_manifest,('sample','truth'))}
        for pred in read_tsv(args.predictions_manifest,('sample','tool')):
            if pred['sample'] not in truth: raise KeyError(f"No truth for {pred['sample']}")
            m=score(truth[pred['sample']],pred); rows.append({'sample':pred['sample'],'tool':pred['tool'],**m})
    else:
        if not args.truth: ap.error('single-sample mode requires --truth')
        pred={'retained_bam':args.retained_bam or '','rrna_r1':args.rrna_fastq or '','non_rrna_r1':args.non_rrna_fastq or ''}
        if not pred['retained_bam'] and not (pred['rrna_r1'] and pred['non_rrna_r1']): ap.error('provide --retained-bam or both FASTQ classification files')
        rows=[{'sample':args.sample,'tool':'custom',**score(args.truth,pred)}]
    p=Path(args.output);p.parent.mkdir(parents=True,exist_ok=True);fields=['sample','tool','TP','FP','TN','FN','sensitivity','specificity','precision','recall','F1','accuracy','false_positive_rate','MCC']
    with p.open('w',newline='') as fh:
        w=csv.DictWriter(fh,fieldnames=fields,delimiter='\t');w.writeheader()
        for r in rows:w.writerow({k:(f'{v:.6f}' if isinstance(v,float) else v) for k,v in r.items()})
    print(p)
if __name__=='__main__':main()
