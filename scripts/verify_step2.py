#!/usr/bin/env python3
"""Validate final matrices against all shards and the independent STEP1 counts."""
import csv,gzip,json
from pathlib import Path
import numpy as np
from step2_pipeline import OUT,S0,S1,MANIFEST,CHROMS,STATS,REGIONS,sha

def verify():
 config=json.loads((OUT/'config.json').read_text());assert sha(config['candidate_source'])==config['candidate_sha256']
 assert sha(S0/'COMPLETE.json')==config['step0_complete_sha256'];assert sha(S1/'COMPLETE.json')==config['step1_complete_sha256']
 for path,old in json.loads((S0/'input_metadata.json').read_text()).items():
  p=Path(path)
  if old is None:assert not p.exists(),path
  else:assert {'size':p.stat().st_size,'mtime_ns':p.stat().st_mtime_ns}==old,path
 codes=np.fromfile(OUT/'candidate_codes.u32',dtype='<u4');n=len(codes);allcodes=np.fromfile(S1/'preliminary_candidates.u32',dtype='<u4');idx=np.searchsorted(allcodes,codes);assert np.array_equal(allcodes[idx],codes)
 A=np.load(OUT/'matrices/matrix_A_counts.npy',mmap_mode='r');B=np.load(OUT/'matrices/matrix_B_counts.npy',mmap_mode='r');AD=np.load(OUT/'matrices/matrix_A_density_per_M.npy',mmap_mode='r');BD=np.load(OUT/'matrices/matrix_B_density_per_M.npy',mmap_mode='r')
 assert A.shape==(n,len(MANIFEST)) and B.shape==(n,len(MANIFEST)*2);assert A.dtype==B.dtype==np.dtype('<u8')
 with np.load(OUT/'matrices/availability_and_denominators.npz') as z:availability=z['available'];den=z['valid_starts']
 pooled=np.zeros((n,len(CHROMS),4),dtype=np.uint64)
 for hi,r in enumerate(MANIFEST):
  h=r['hap']
  with np.load(OUT/'counts'/(h+'.npz')) as z:oriented=z['counts'];chrs=z['chromosomes'].tolist();d=z['valid_starts'];bp=z['bp'];assert z['regions'].tolist()==REGIONS
  c=oriented.sum(axis=3);total=c.sum(axis=1);ci=[CHROMS.index(ch) for ch in chrs]
  assert np.array_equal(np.flatnonzero(availability[hi]),ci);assert np.array_equal(den[hi,ci],d)
  for j,ch in enumerate(chrs):
   for ri,f in enumerate(REGIONS):assert int(d[j,ri])==int(STATS[(h,ch,f)]['valid_16mer_starts']) and int(bp[j,ri])==int(STATS[(h,ch,f)]['bp'])
  assert np.array_equal(A[:,hi],total[:,2]);assert np.array_equal(B[:,hi*2:hi*2+2],total[:,:2])
  with np.errstate(divide='ignore',invalid='ignore'):
   assert np.allclose(AD[:,hi],total[:,2]/d[:,2].sum()*1e6,equal_nan=True,rtol=1e-12,atol=0)
   assert np.allclose(BD[:,hi*2:hi*2+2],total[:,:2]/d[:,:2].sum(axis=0)*1e6,equal_nan=True,rtol=1e-12,atol=0)
  with np.load(S1/'counts'/(h+'.npz')) as z:old=z['counts'][idx,:3,:]
  assert np.array_equal(oriented[:,:,:3,:].sum(axis=1),old),h
  assert np.all(c[:,:,2]>=c[:,:,0]+c[:,:,1]+c[:,:,3])
  residual=c[:,:,2]-c[:,:,0]-c[:,:,1]-c[:,:,3]
  assert int(residual.sum())<=int(d[:,2].sum())-int(d[:,0].sum())-int(d[:,1].sum())-int(d[:,3].sum())
  if r['sample']!='CHM13':pooled[:,ci,:]+=c
  if (hi+1)%50==0:print('verify',hi+1,len(MANIFEST),flush=True)
 with np.load(OUT/'matrices/pooled_chromosome_counts.npz') as z:assert np.array_equal(z['counts'],pooled)
 # Independently check the readable candidate summary against pooled matrices.
 totals=pooled.sum(axis=1);reference=np.array([r['sample']=='CHM13' for r in MANIFEST]);prev=(A[:,~reference]>0).sum(axis=1)
 with gzip.open(OUT/'candidate_metrics.tsv.gz','rt') as f:
  row_count=0
  for i,row in enumerate(csv.DictReader(f,delimiter='\t')):
   assert row['kmer_id']==f'k16_{codes[i]:08x}'
   for j,key in enumerate(['HSat2_count','HSat3_count','HSat23_union_count']):assert int(row[key])==int(totals[i,j])
   assert int(row['union_prevalence_haps'])==int(prev[i])
   denom=int(totals[i,0])+int(totals[i,1]);expected=int(totals[i,0])/denom if denom else np.nan
   assert np.isclose(float(row['HSat2_count_fraction']),expected,equal_nan=True,rtol=1e-12,atol=0)
   row_count+=1
  assert row_count==n
 summary=json.loads((OUT/'summary.json').read_text());assert summary['observed_hap_chromosome_pairs']==int(availability.sum())==8531
 assert summary['candidate_count']==n==config['candidate_count']
 assert len(list((OUT/'plots').glob('*.png')))==3
 result={'status':'PASS','candidate_count':n,'sequence_sets':len(MANIFEST),'source_metadata_unchanged':True,'all_shards_match_step1_both_orientations':True,'all_matrices_match_shards':True,'density_normalization_verified':True,'missing_chromosome_mask_verified':True,'pooled_counts_exclude_CHM13':True,'candidate_metrics_counts_prevalence_family_fraction_verified':True,'family_junction_accounting_verified':True}
 (OUT/'validation_summary.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2),flush=True)
if __name__=='__main__':verify()
