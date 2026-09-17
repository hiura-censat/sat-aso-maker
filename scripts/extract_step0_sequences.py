#!/usr/bin/env python3
"""Stream one STEP 0 region view as FASTA; pipe to gzip if desired."""
import argparse, csv, sys
from collections import defaultdict
from step0_preprocess import OUT, records
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--hap',required=True)
p.add_argument('--region',required=True,choices=['HSat2','HSat3','HSat23','ambiguous','background'])
a=p.parse_args()
with (OUT/'manifest.tsv').open() as f:
    match=[r for r in csv.DictReader(f,delimiter='\t') if r['hap']==a.hap and r['status']=='PASS']
if len(match)!=1: p.error('hap is not eligible in manifest.tsv')
regions=defaultdict(list)
with (OUT/'regions'/a.hap/(a.region+'.bed')).open() as f:
    for line in f:
        c,s,e,*_=line.split(); regions[c].append((int(s),int(e)))
for c,seq in records(match[0]['fasta']):
    for s,e in regions[c]:
        sys.stdout.buffer.write(f'>{a.hap}|{c}:{s}-{e}|{a.region}|strand=+\n'.encode())
        for offset in range(s,e,80): sys.stdout.buffer.write(seq[offset:min(offset+80,e)]+b'\n')
