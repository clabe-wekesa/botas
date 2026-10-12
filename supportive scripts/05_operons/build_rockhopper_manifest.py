#!/usr/bin/env python3
import argparse,csv,sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'lib'))
from manifest import read_tsv

def main():
    ap=argparse.ArgumentParser(description='Join samples.tsv to prepared Rockhopper reference directories.')
    ap.add_argument('--samples',required=True,help='TSV columns: sample,r1,r2')
    ap.add_argument('--reference-dir',required=True,help='Root created by make_rockhopper_annotations.py')
    ap.add_argument('--output',required=True);args=ap.parse_args()
    root=Path(args.reference_dir);rows=[]
    for x in read_tsv(args.samples,('sample','r1','r2')):
        d=root/x['sample']
        if not d.is_dir(): raise FileNotFoundError(d)
        rows.append(dict(sample=x['sample'],r1=str(Path(x['r1']).resolve()),r2=str(Path(x['r2']).resolve()),rockhopper_reference=str(d.resolve())))
    p=Path(args.output);p.parent.mkdir(parents=True,exist_ok=True)
    with p.open('w',newline='') as fh:w=csv.DictWriter(fh,fieldnames=['sample','r1','r2','rockhopper_reference'],delimiter='\t');w.writeheader();w.writerows(rows)
    print(p)
if __name__=='__main__':main()
