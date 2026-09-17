#!/usr/bin/env python3
"""Export selected STEP 1 candidates for one hap as a tab-separated table."""
import argparse,csv,sys
import numpy as np
from step1_pipeline import OUT,REGIONS,sequence
p=argparse.ArgumentParser(description=__doc__);p.add_argument('--hap',required=True);p.add_argument('--kmer',help='Optional 16-mer; either orientation is accepted');a=p.parse_args()
codes=np.fromfile(OUT/'preliminary_candidates.u32',dtype='<u4')
if a.kmer:
 seq=a.kmer.upper()
 if len(seq)!=16 or set(seq)-set('ACGT'):p.error('--kmer must contain exactly 16 ACGT bases')
 seq=min(seq,seq.translate(str.maketrans('ACGT','TGCA'))[::-1]);code=0
 for c in seq:code=(code<<2)|'ACGT'.index(c)
 idx=np.searchsorted(codes,code)
 if idx>=len(codes) or codes[idx]!=code:p.error('Not in the preliminary candidate set; background was not counted for this k-mer.')
 indices=[idx]
else:
 chosen=np.fromfile(OUT/'candidate_codes.u32',dtype='<u4');indices=np.searchsorted(codes,chosen)
file=OUT/'counts'/(a.hap+'.npz')
if not file.exists():p.error('Unknown or incomplete hap')
with np.load(file) as z:counts=z['counts']
with np.load(OUT/'enrichment_by_hap'/(a.hap+'.npz')) as z:t=z['target_density'];b=z['background_density'];e=z['log2_enrichment']
w=csv.writer(sys.stdout,delimiter='\t');w.writerow(['hap','kmer_id','canonical_kmer']+[f'{f}_{o}' for f in REGIONS for o in ['forward','reverse']]+['target_density','background_density','log2_enrichment'])
for i in indices:w.writerow([a.hap,f'k16_{int(codes[i]):08x}',sequence(codes[i])]+[int(x) for x in counts[i].ravel()]+[t[i],b[i],e[i]])
