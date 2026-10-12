#!/usr/bin/env python3
import argparse,csv,re
POSITIVE=[r'northern blot',r'rt[\s-]?pcr',r'rnase protection',r's1 nuclease',r's1 mapping',r'transcript(?:ion)?\s+analysis',r'transcript(?:ion)?\s+unit',r'polycistronic',r'transcriptional fusion']
NEGATIVE=[r'opposite direction',r'opposite strand',r'separate transcript',r'separately transcribed',r'independent transcript',r'independently transcribed',r'genome analysis',r'prediction']
def strict(text):
    e=text.lower(); return not any(re.search(p,e) for p in NEGATIVE) and any(re.search(p,e) for p in POSITIVE)
def main():
    ap=argparse.ArgumentParser(description='Filter DBTBS operons to strong experimental evidence and emit adjacent truth pairs.')
    ap.add_argument('--input',required=True);ap.add_argument('--operons-out',required=True);ap.add_argument('--pairs-out',required=True);args=ap.parse_args()
    with open(args.input,newline='') as fh: rows=[r for r in csv.DictReader(fh,delimiter='\t') if len([x for x in r['genes'].split(',') if x.strip()])>=2 and strict(r['evidence'])]
    if not rows: raise RuntimeError('No strict operons retained')
    with open(args.operons_out,'w',newline='') as fh: w=csv.DictWriter(fh,fieldnames=rows[0].keys(),delimiter='\t');w.writeheader();w.writerows(rows)
    with open(args.pairs_out,'w',newline='') as fh:
        w=csv.writer(fh,delimiter='\t');w.writerow(['gene1','gene2'])
        seen=set()
        for r in rows:
            g=[x.strip() for x in r['genes'].split(',') if x.strip()]
            for a,b in zip(g,g[1:]):
                if (a,b) not in seen: w.writerow([a,b]);seen.add((a,b))
if __name__=='__main__':main()
