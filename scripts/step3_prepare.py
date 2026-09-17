#!/usr/bin/env python3
import csv,json,hashlib,time
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'step3';S2=ROOT/'step2'
def tsv(path,head,rows):
 with path.open('w') as f:
  w=csv.writer(f,delimiter='\t');w.writerow(head);w.writerows(rows)
def prepare():
 for sub in ['data','clustering','plots','candidates','mismatch','logs']:(OUT/sub).mkdir(parents=True,exist_ok=True)
 haps=list(csv.DictReader((S2/'haplotypes.tsv').open(),delimiter='\t'));codes=np.fromfile(S2/'candidate_codes.u32',dtype='<u4');n=len(codes);nh=len(haps)
 with np.load(S2/'matrices/availability_and_denominators.npz') as z:den=z['valid_starts'];avail=z['available'];ref=z['is_reference'];ch=z['chromosomes']
 core=[np.flatnonzero((den[~ref,:,f]>0).mean(axis=0)>=.85) for f in range(2)]
 config={'seed':42,'core_chromosome_min_measurable_fraction':.85,'core_chromosomes':{f:ch[core[i]].tolist() for i,f in enumerate(['HSat2','HSat3'])},'candidate_count':n,'nonreference_haps':int((~ref).sum()),'chromosomes':ch.tolist(),'source_complete_sha256':hashlib.sha256((S2/'COMPLETE.json').read_bytes()).hexdigest()}
 (OUT/'config.json').write_text(json.dumps(config,indent=2)+'\n')
 tsv(OUT/'coverage.tsv',['chromosome','available_haps','HSat2_measurable_haps','HSat3_measurable_haps','HSat2_core','HSat3_core'],((str(c),int(avail[~ref,i].sum()),int((den[~ref,i,0]>0).sum()),int((den[~ref,i,1]>0).sum()),int(i in core[0]),int(i in core[1])) for i,c in enumerate(ch)))
 if (OUT/'data/PREPARED.json').exists():print('PREPARED cached',flush=True);return
 B=np.load(S2/'matrices/matrix_B_counts.npy',mmap_mode='r');BD=np.load(S2/'matrices/matrix_B_density_per_M.npy',mmap_mode='r')
 pools=[]
 for f in range(2):
  c=np.asarray(B[:,f::2][:,~ref]);x=np.log1p(np.asarray(BD[:,f::2][:,~ref]));valid=(c>=10).sum(axis=1)>=29;v=np.nanvar(x,axis=1);v[~valid]=-1
  pools.append(np.argsort(-v,kind='stable')[:5000]);del c,x,v
 pool=np.array(pools);np.save(OUT/'data/chromosome_feature_indices.npy',pool)
 C=np.lib.format.open_memmap(OUT/'data/core_counts.npy',mode='w+',dtype='<u8',shape=(nh,2,n))
 CC=np.lib.format.open_memmap(OUT/'data/chromosome_feature_counts.npy',mode='w+',dtype='<u4',shape=(nh,len(ch),2,5000));CC[:]=0
 coreden=np.zeros((nh,2),dtype=np.uint64);complete=np.zeros((nh,2),dtype=bool)
 for hi,h in enumerate(haps):
  with np.load(S2/'counts'/(h['hap']+'.npz')) as z:c=z['counts'].sum(axis=3);chs=z['chromosomes'].tolist()
  ci=np.array([list(ch).index(s) for s in chs])
  for f in range(2):
   use=np.isin(ci,core[f]);C[hi,f]=c[:,use,f].sum(axis=1);coreden[hi,f]=den[hi,core[f],f].sum();complete[hi,f]=np.all(den[hi,core[f],f]>0)
   selected=c[pool[f],:,f].T;assert selected.max()<2**32;CC[hi,ci,f,:]=selected
  if (hi+1)%25==0:print('prepare',hi+1,nh,flush=True)
 C.flush();CC.flush();np.savez_compressed(OUT/'data/metadata.npz',codes=codes,haplotypes=np.array([h['hap'] for h in haps]),is_reference=ref,chromosomes=ch,valid_starts=den,available=avail,core_denominators=coreden,core_complete=complete)
 (OUT/'data/PREPARED.json').write_text(json.dumps({'status':'PASS','core_complete_nonreference':complete[~ref].sum(axis=0).tolist(),'candidate_count':n},indent=2)+'\n')
 print('PREPARED',complete[~ref].sum(axis=0),flush=True)
if __name__=='__main__':prepare()
