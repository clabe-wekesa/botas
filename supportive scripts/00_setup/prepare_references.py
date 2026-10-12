#!/usr/bin/env python3
"""Download and validate version-pinned RefSeq genome/GFF3 pairs from a TSV manifest.
Requires Python >=3.9 and curl. No pip, conda or NCBI account required.
Usage: python3 prepare_references.py --manifest references.tsv --out references
Includes every replicon supplied in each assembly. Does not rotate references.
Use these exact FASTA files for all aligners. Supply your previous simulation
references separately if sequence versions/origins differ.
"""
import argparse
import csv
import datetime
import gzip
import hashlib
import re
import shutil
import subprocess
import sys
from pathlib import Path
from urllib.parse import unquote



def fetch(url, path):
 path.parent.mkdir(parents=True, exist_ok=True)
 subprocess.run(['curl','--fail','--location','--retry','5','--retry-delay','3',
                 '--connect-timeout','30','--silent','--show-error',
                 '--output',str(path),url],check=True)

def digest(path, algorithm='sha256'):
 h=hashlib.new(algorithm)
 with path.open('rb') as f:
  for b in iter(lambda:f.read(1024*1024),b''): h.update(b)
 return h.hexdigest()

def checked_download(base, name, md5s, directory):
 p=directory/name
 if name not in md5s: raise ValueError('Missing NCBI MD5 entry: '+name)
 if p.exists() and digest(p,'md5')==md5s[name]:
  print('  Verified existing '+name,flush=True)
  return p
 part=p.with_suffix(p.suffix+'.part')
 fetch(base+'/'+name,part)
 if digest(part,'md5')!=md5s[name]: raise ValueError('MD5 mismatch: '+str(part))
 with gzip.open(part,'rb') as f:
  while f.read(1024*1024): pass
 part.replace(p)
 return p

def fasta_records(path):
 name=None; chunks=[]
 with gzip.open(path,'rt') as f:
  for line in f:
   if line.startswith('>'):
    if name is not None: yield name,''.join(chunks)
    name=line[1:].split()[0];chunks=[]
   elif line.strip():
    if name is None: raise ValueError('Sequence before FASTA header')
    chunks.append(line.strip())
 if name is not None: yield name,''.join(chunks)

def prepare_fasta(source, target):
 lengths={}; rows=[]
 part=target.with_suffix('.fna.tmp')
 with part.open('wb') as f:
  for name,sequence in fasta_records(source):
   if name in lengths or not sequence: raise ValueError('Duplicate/empty FASTA sequence: '+name)
   if not re.fullmatch('[ACGTRYSWKMBDHVNacgtryswkmbdhvn]+',sequence):
    raise ValueError('Unexpected DNA character in '+name)
   lengths[name]=len(sequence)
   f.write(('>'+name+'\n').encode());offset=f.tell()
   for i in range(0,len(sequence),60): f.write((sequence[i:i+60]+'\n').encode())
   width=min(60,len(sequence));rows.append((name,len(sequence),offset,width,width+1))
 if not lengths: raise ValueError('Empty FASTA')
 part.replace(target)
 index=Path(str(target)+'.fai')
 index.write_text(''.join('\t'.join(map(str,r))+'\n' for r in rows))
 return lengths

def prepare_gff(source, target, lengths):
 counts={}; regions={}; seen=set();has_version=False
 part=target.with_suffix('.gff3.tmp')
 with gzip.open(source,'rt') as src, part.open('w') as dst:
  for line in src:
   if line.startswith('##FASTA'): break
   if line.startswith('##gff-version 3'): has_version=True
   if line.startswith('##sequence-region'):
    fields=line.split()
    if len(fields)!=4 or fields[1] not in lengths: raise ValueError('Invalid GFF sequence-region')
    if int(fields[2])!=1 or int(fields[3])!=lengths[fields[1]]: raise ValueError('GFF sequence-region length mismatch')
   if line.startswith('#') or not line.strip(): dst.write(line);continue
   fields=line.rstrip('\n').split('\t')
   if len(fields)!=9: raise ValueError('GFF must have nine columns')
   seq,kind=fields[0],fields[2];start,end=int(fields[3]),int(fields[4])
   if seq not in lengths: raise ValueError('GFF sequence ID absent from FASTA: '+seq)
   if not 1<=start<=end<=lengths[seq]: raise ValueError('GFF coordinates out of bounds: '+seq)
   seen.add(seq);counts[kind]=counts.get(kind,0)+1
   if kind=='region':
    attrs=dict((k,unquote(v)) for k,v in (a.split('=',1) for a in fields[8].split(';') if '=' in a))
    regions[seq]=attrs
   dst.write(line)
 if not has_version: raise ValueError('Not declared GFF3')
 if not counts.get('gene'): raise ValueError('Annotation has no gene features')
 missing=set(lengths)-seen
 if missing: raise ValueError('Replicons have no annotation: '+','.join(sorted(missing)))
 part.replace(target)
 return counts,regions

def write_tsv(path,rows,columns):
 part=path.with_suffix('.tmp')
 with part.open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=columns,delimiter='\t');w.writeheader();w.writerows(rows)
 part.replace(path)

def main():
 ap=argparse.ArgumentParser(description=__doc__,formatter_class=argparse.RawDescriptionHelpFormatter)
 ap.add_argument('--manifest',required=True,help='TSV columns: sample, assembly, assembly_name')
 ap.add_argument('--out',required=True)
 args=ap.parse_args()
 import csv
 with open(args.manifest,newline='') as fh:
  rr=csv.DictReader(fh,delimiter='\t')
  required={'sample','assembly','assembly_name'}
  if not required.issubset(rr.fieldnames or []): ap.error('reference manifest must contain sample, assembly, assembly_name')
  references=[(r['sample'].strip(),r['assembly'].strip(),r['assembly_name'].strip()) for r in rr if r['sample'].strip()]
 if not shutil.which('curl'): ap.error('curl is required: sudo apt install curl')
 root=Path(args.out).resolve();root.mkdir(parents=True,exist_ok=True)
 manifest=[];replicons=[];errors=[]
 for species,accession,assembly in references:
  print('\n'+species+' '+accession,flush=True)
  try:
   digits=accession.split('_')[1].split('.')[0]
   prefix=accession.split('_')[0]
   folder=accession+'_'+assembly
   base=f'https://ftp.ncbi.nlm.nih.gov/genomes/all/{prefix}/{digits[:3]}/{digits[3:6]}/{digits[6:9]}/{folder}'
   directory=root/species;raw=directory/'original';raw.mkdir(parents=True,exist_ok=True)
   md5path=raw/'md5checksums.txt';fetch(base+'/md5checksums.txt',md5path)
   md5s={}
   for line in md5path.read_text().splitlines():
    values=line.split()
    if len(values)==2:md5s[Path(values[1].lstrip('*')).name]=values[0]
   fasta=checked_download(base,folder+'_genomic.fna.gz',md5s,raw)
   gff=checked_download(base,folder+'_genomic.gff.gz',md5s,raw)
   prepared_fasta=directory/(species+'.fna');prepared_gff=directory/(species+'.gff3')
   lengths=prepare_fasta(fasta,prepared_fasta)
   counts,regions=prepare_gff(gff,prepared_gff,lengths)
   for seq,length in lengths.items():
    attrs=regions.get(seq,{})
    replicons.append(dict(species=species,assembly=accession,seqid=seq,length=length,
                         circular=attrs.get('Is_circular','unknown'),
                         chromosome=attrs.get('chromosome',''),plasmid=attrs.get('plasmid','')))
   manifest.append(dict(species=species,assembly=accession,fasta=str(prepared_fasta),gff3=str(prepared_gff),
                        fai=str(prepared_fasta)+'.fai',replicons=len(lengths),bases=sum(lengths.values()),
                        genes=counts.get('gene',0),CDS=counts.get('CDS',0),
                        fasta_sha256=digest(prepared_fasta),gff3_sha256=digest(prepared_gff),
                        source_url=base,retrieved_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()))
   print(f'  READY: {len(lengths)} replicons; {sum(lengths.values()):,} bp; {counts.get("gene",0):,} genes',flush=True)
  except Exception as exc:
   errors.append(species+': '+str(exc));print('  FAILED: '+str(exc),file=sys.stderr,flush=True)
 write_tsv(root/'references.tsv',manifest,['species','assembly','fasta','gff3','fai','replicons','bases','genes','CDS','fasta_sha256','gff3_sha256','source_url','retrieved_utc'])
 write_tsv(root/'replicons.tsv',replicons,['species','assembly','seqid','length','circular','chromosome','plasmid'])
 (root/'errors.log').write_text('\n'.join(errors))
 print(f'\nPrepared {len(manifest)} references. Manifest: {root / "references.tsv"}')
 print('All included replicons retained; sequence IDs and coordinates preserved. FASTA .fai indexes created.')
 print('Compare these reference sequences with your simulation files before claiming identical references.')
 return 1 if errors else 0

if __name__=='__main__':
 try:sys.exit(main())
 except KeyboardInterrupt:
  print('\nStopped. Rerun to reuse checksum-verified downloads.',file=sys.stderr);sys.exit(130)

