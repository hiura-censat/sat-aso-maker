#!/usr/bin/env python3
import csv,json,hashlib,time
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'step3';S2=ROOT/'step2'
def tsv(path,head,rows):
 with path.open('w') as f:
  w=csv.writer(f,delimiter='\t');w.writerow(head);w.writerows(rows)
def choose_feature_pool(pooled_counts,width=5000):
 n,nchr,_=pooled_counts.shape;width=min(width,n);pool=np.empty((nchr,2,width),dtype=np.int32)
 for ci in range(nchr):
  for f in range(2):
   observed=pooled_counts[:,ci,f]
   ordered=np.argsort(-observed.astype(np.int64),kind='stable');positive=ordered[observed[ordered]>0]
   common=positive[:width//2];remainder=positive[width//2:]
   diverse=remainder[np.linspace(0,len(remainder)-1,min(width-len(common),len(remainder)),dtype=int)] if len(remainder) else np.array([],dtype=int)
   chosen=np.concatenate([common,diverse]);unused=np.ones(n,dtype=bool);unused[chosen]=False
   pool[ci,f]=np.concatenate([chosen,ordered[unused[ordered]][:width-len(chosen)]])
 return pool
def prepare():
 for sub in ['data','clustering','plots','candidates','mismatch','logs']:(OUT/sub).mkdir(parents=True,exist_ok=True)
 haps=list(csv.DictReader((S2/'haplotypes.tsv').open(),delimiter='\t'));codes=np.fromfile(S2/'candidate_codes.u32',dtype='<u4');n=len(codes);nh=len(haps)
 with np.load(S2/'matrices/availability_and_denominators.npz') as z:den=z['valid_starts'];avail=z['available'];ref=z['is_reference'];ch=z['chromosomes']
 config={'seed':42,'clustering_unit':'family_by_chromosome','minimum_measurable_haps':100,'candidate_count':n,'nonreference_haps':int((~ref).sum()),'chromosomes':ch.tolist(),'source_complete_sha256':hashlib.sha256((S2/'COMPLETE.json').read_bytes()).hexdigest()}
 (OUT/'config.json').write_text(json.dumps(config,indent=2)+'\n')
 tsv(OUT/'coverage.tsv',['chromosome','available_haps','HSat2_measurable_haps','HSat3_measurable_haps'],((str(c),int(avail[~ref,i].sum()),int((den[~ref,i,0]>0).sum()),int((den[~ref,i,1]>0).sum())) for i,c in enumerate(ch)))
 if (OUT/'data/PREPARED.json').exists():print('PREPARED cached',flush=True);return
 with np.load(S2/'matrices/pooled_chromosome_counts.npz') as z:
  assert np.array_equal(z['codes'],codes) and np.array_equal(z['chromosomes'],ch)
  pool=choose_feature_pool(z['counts']);width=pool.shape[2]
 np.save(OUT/'data/chromosome_feature_indices.npy',pool)
 CC=np.lib.format.open_memmap(OUT/'data/chromosome_feature_counts.npy',mode='w+',dtype='<u4',shape=(nh,len(ch),2,width));CC[:]=0
 for hi,h in enumerate(haps):
  with np.load(S2/'counts'/(h['hap']+'.npz')) as z:c=z['counts'].sum(axis=3);chs=z['chromosomes'].tolist()
  ci=np.array([list(ch).index(s) for s in chs])
  for local,global_ci in enumerate(ci):
   for f in range(2):
    selected=c[pool[global_ci,f],local,f];assert selected.max()<2**32;CC[hi,global_ci,f,:]=selected
  if (hi+1)%25==0:print('prepare',hi+1,nh,flush=True)
 CC.flush();np.savez_compressed(OUT/'data/metadata.npz',codes=codes,haplotypes=np.array([h['hap'] for h in haps]),is_reference=ref,chromosomes=ch,valid_starts=den,available=avail)
 (OUT/'data/PREPARED.json').write_text(json.dumps({'status':'PASS','measurable_nonreference_by_chromosome_family':(den[~ref,:,:2]>0).sum(axis=0).tolist(),'candidate_count':n},indent=2)+'\n')
 print('PREPARED',flush=True)
if __name__=='__main__':prepare()
