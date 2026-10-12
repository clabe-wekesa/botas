#!/usr/bin/env python3
from __future__ import annotations
import csv
from pathlib import Path

def read_tsv(path, required=()):
    p=Path(path)
    if not p.is_file(): raise FileNotFoundError(p)
    with p.open(newline='') as fh:
        r=csv.DictReader(fh, delimiter='\t')
        fields=r.fieldnames or []
        missing=[x for x in required if x not in fields]
        if missing: raise ValueError(f"{p}: missing columns: {', '.join(missing)}")
        rows=[]
        for n,row in enumerate(r,2):
            if not any((v or '').strip() for v in row.values()): continue
            row={k:(v or '').strip() for k,v in row.items()}
            for key in required:
                if not row[key]: raise ValueError(f"{p}:{n}: empty required field {key}")
            rows.append(row)
    return rows

def ensure_file(path,label='file'):
    p=Path(path).expanduser()
    if not p.is_file(): raise FileNotFoundError(f"{label} not found: {p}")
    return p
