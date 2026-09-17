#!/usr/bin/env python3
import csv,itertools,json
import numpy as np
from step3_prepare import OUT,tsv
def revcode(code):
 out=0
 for _ in range(16):out=(out<<2)|((code&3)^3);code>>=2
 return out
def seq(code):return ''.join('ACGT'[(int(code)>>s)&3] for s in range(30,-1,-2))
def encode(s):
 code=0
 for b in s:code=(code<<2)|'ACGT'.index(b)
 return code
def neighbors(code):
 variants={min(code,revcode(code))}
 for i in range(16):
  for b in range(1,4):
   v=code^(b<<(2*i));variants.add(min(v,revcode(v)))
 for i,j in itertools.combinations(range(16),2):
  for b in range(1,4):
   for d in range(1,4):
    v=code^(b<<(2*i))^(d<<(2*j));variants.add(min(v,revcode(v)))
 return variants
def distance_masks(query,variant):
 q=seq(query);v=seq(variant);r=seq(revcode(variant));masks=[sum(1<<i for i,(a,b) in enumerate(zip(q,s)) if a!=b) for s in [v,r]];dist=min(m.bit_count() for m in masks);best=sorted(set(m for m in masks if m.bit_count()==dist));return dist,best[0],best[-1]
def build():
 queries=list(csv.DictReader((OUT/'candidates/mismatch_queries.tsv').open(),delimiter='\t'));qcodes=[encode(r['canonical_kmer']) for r in queries];mapping=[]
 for qi,code in enumerate(qcodes):
  for variant in sorted(neighbors(code)):
   d,m1,m2=distance_masks(code,variant);assert d<=2;mapping.append((variant,qi,d,m1,m2))
 codes=np.array(sorted({v[0] for v in mapping}),dtype='<u4');codes.tofile(OUT/'mismatch/variant_codes.u32');lookup={int(c):i for i,c in enumerate(codes)}
 records=np.array([(lookup[v],q,d,m1,m2) for v,q,d,m1,m2 in mapping],dtype=np.uint32)
 np.savez_compressed(OUT/'mismatch/variant_mapping.npz',records=records,query_codes=np.array(qcodes,dtype='<u4'),variant_codes=codes)
 tsv(OUT/'mismatch/variant_catalog.tsv',['variant_index','variant_canonical_kmer','query_index','query_kmer_id','minimum_Hamming_distance','mismatch_positions_1based','alternative_tied_positions_1based'],((lookup[v],seq(v),q,queries[q]['kmer_id'],d,','.join(str(i+1) for i in range(16) if m1&(1<<i)),','.join(str(i+1) for i in range(16) if m2&(1<<i)) if m1!=m2 else '') for v,q,d,m1,m2 in mapping))
 (OUT/'mismatch/config.json').write_text(json.dumps({'queries':len(queries),'distinct_canonical_variants':len(codes),'query_variant_relations':len(records),'max_substitutions':2,'indels':False,'distance':'minimum over both orientations; once per genomic window per query','position_ties':'equal weight to distinct tied orientation masks'},indent=2)+'\n');print('VARIANTS',len(queries),len(codes),flush=True)
if __name__=='__main__':build()
