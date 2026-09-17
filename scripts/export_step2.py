#!/usr/bin/env python3
"""Export one candidate across all chromosome/family/hap columns, preserving NA."""
import argparse,csv,sys
import numpy as np
from step2_pipeline import OUT,MANIFEST,CHROMS,REGIONS

def export(kmer,hap=None,stream=sys.stdout):
 kmer=kmer.upper()
 if len(kmer)!=16 or set(kmer)-set('ACGT'):raise ValueError('kmer must be 16 A/C/G/T bases')
 rev=kmer.translate(str.maketrans('ACGT','TGCA'))[::-1];canonical=min(kmer,rev);code=0
 for b in canonical:code=(code<<2)|'ACGT'.index(b)
 codes=np.fromfile(OUT/'candidate_codes.u32',dtype='<u4');i=int(np.searchsorted(codes,code))
 if i==len(codes) or codes[i]!=code:raise ValueError('kmer is not in the STEP2 candidate set; not a zero count')
 jobs=[r for r in MANIFEST if hap is None or r['hap']==hap]
 if not jobs:raise ValueError('unknown hap')
 w=csv.writer(stream,delimiter='\t');w.writerow(['kmer_id','canonical_kmer','hap','is_reference','chromosome','region','available','canonical_forward_count','reverse_complement_count','total_count','valid_16mer_starts','density_per_M'])
 for r in jobs:
  with np.load(OUT/'counts'/(r['hap']+'.npz')) as z:counts=z['counts'][i];chromosomes=z['chromosomes'].tolist();den=z['valid_starts']
  for ch in CHROMS:
   for ri,region in enumerate(REGIONS):
    base=[f'k16_{code:08x}',canonical,r['hap'],int(r['sample']=='CHM13'),ch,region]
    if ch not in chromosomes:w.writerow(base+[0]+['NA']*5);continue
    ci=chromosomes.index(ch);f,b=map(int,counts[ci,ri]);d=int(den[ci,ri]);w.writerow(base+[1,f,b,f+b,d,(f+b)/d*1e6 if d else 'NA'])
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--kmer',required=True);p.add_argument('--hap');a=p.parse_args();export(a.kmer,a.hap)
