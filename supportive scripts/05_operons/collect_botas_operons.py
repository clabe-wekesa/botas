#!/usr/bin/env python3
import argparse,csv,sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'lib'))
from manifest import read_tsv,ensure_file

def main():
    ap=argparse.ArgumentParser(description='Combine BOTAS operon TSVs listed in a manifest.')
    ap.add_argument('--manifest',required=True,help='TSV columns: sample,operons');ap.add_argument('--output',required=True);args=ap.parse_args()
    rows=read_tsv(args.manifest,('sample','operons')); allrows=[];fields=None
    for x in rows:
        with open(ensure_file(x['operons'],'operon TSV'),newline='') as fh:
            rr=csv.DictReader(fh,delimiter='\t');fields=fields or ['sample',*(rr.fieldnames or [])]
            allrows.extend({'sample':x['sample'],**r} for r in rr)
    p=Path(args.output);p.parent.mkdir(parents=True,exist_ok=True)
    with p.open('w',newline='') as fh:w=csv.DictWriter(fh,fieldnames=fields,delimiter='\t');w.writeheader();w.writerows(allrows)
    print(p)
if __name__=='__main__':main()
