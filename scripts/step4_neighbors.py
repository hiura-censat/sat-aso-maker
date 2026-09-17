#!/usr/bin/env python3
import csv,json
import numpy as np
from step3_neighbors import encode,neighbors,distance_masks,seq
from step3_prepare import ROOT
OUT=ROOT/'step4';M=OUT/'mismatch'
def main():
 rows=list(csv.DictReader((M/'queries.tsv').open(),delimiter='\t'));new=[r for r in rows if r['query_origin']=='STEP4'];qcodes=[encode(r['canonical_kmer']) for r in new];mapping=[]
 for qi,code in enumerate(qcodes):
  for variant in sorted(neighbors(code)):
   d,m1,m2=distance_masks(code,variant);mapping.append((variant,qi,d,m1,m2))
 codes=np.array(sorted({r[0] for r in mapping}),dtype='<u4');codes.tofile(M/'variant_codes.u32');lookup={int(c):i for i,c in enumerate(codes)};records=np.array([(lookup[v],q,d,m1,m2) for v,q,d,m1,m2 in mapping],dtype=np.uint32)
 np.savez_compressed(M/'variant_mapping.npz',records=records,query_codes=np.array(qcodes,dtype='<u4'),variant_codes=codes)
 with (M/'new_queries.tsv').open('w') as f:
  fields=list(rows[0]) if rows else []
  w=csv.DictWriter(f,fieldnames=fields,delimiter='\t');w.writeheader();w.writerows(new)
 (M/'config.json').write_text(json.dumps({'new_queries':len(new),'existing_step3_queries':len(rows)-len(new),'total_queries':len(rows),'distinct_new_canonical_variants':len(codes),'new_query_variant_relations':len(records),'max_substitutions':2,'indels':False},indent=2)+'\n')
 print('STEP4_VARIANTS',len(new),len(codes),flush=True)
if __name__=='__main__':main()
