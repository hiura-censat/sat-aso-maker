#!/usr/bin/env python3
"""Sweep strict E > thresholds using existing exact counts; retain original STEP1."""
import csv, json, math
import numpy as np
from step1_pipeline import OUT, MANIFEST, STATS, enrichment
from workflow_settings import settings

DEST=OUT/'threshold_sweep';DEST.mkdir(exist_ok=True)
SELECT=settings()['selection']
THRESHOLDS=SELECT['sweep_E']
SELECTED=SELECT['selected_E']
POOLED_MIN=SELECT['pooled_min_count']
HAP_MIN=SELECT['hap_min_count']
BROAD_FRACTION=SELECT['broad_hap_fraction']
with np.load(OUT/'pooled_counts_and_scores.npz') as z:
 codes=z['codes'];e=z['log2_enrichment'];counts=z['counts']
 t=counts[:,2,:].sum(axis=1);b=counts[:,3,:].sum(axis=1)
 prev=z['prevalence_haps'];td=z['target_density'];bd=z['background_density']
 del counts
n=len(codes);nh=sum(r['sample']!='CHM13' for r in MANIFEST)
cache=DEST/'local_threshold_counts.npz'
if cache.exists():
 with np.load(cache) as z:
  assert np.array_equal(z['codes'],codes) and np.array_equal(z['thresholds'],THRESHOLDS)
  local=z['qualifying_haps']
else:
 local=np.zeros((len(THRESHOLDS),n),dtype=np.uint16)
 for i,r in enumerate(MANIFEST,1):
  if r['sample']=='CHM13':continue
  h=r['hap']
  with np.load(OUT/'counts'/(h+'.npz')) as z:c=z['counts']
  ht=c[:,2,:].sum(axis=1);hb=c[:,3,:].sum(axis=1)
  _,_,he=enrichment(ht,hb,STATS[(h,'HSat23')],STATS[(h,'background')])
  for j,threshold in enumerate(THRESHOLDS):local[j]+=((ht>=HAP_MIN)&(he>threshold)).astype(np.uint16)
  if i%25==0:print('sweep',i,len(MANIFEST),flush=True)
 np.savez_compressed(cache,codes=codes,thresholds=THRESHOLDS,qualifying_haps=local)
assert np.all(local[1:]<=local[:-1])
rows=[]
for j,threshold in enumerate(THRESHOLDS):
 pooled=(t>=POOLED_MIN)&(e>threshold);hap=local[j]>0
 for route,mask in [('pooled',pooled),('hap_only',hap&~pooled),('union',pooled|hap)]:
  k=int(mask.sum())
  rows.append([threshold,2**threshold,route,k,int((mask&(b==0)).sum()),int((mask&(prev>=math.ceil(BROAD_FRACTION*nh))).sum()),int((mask&(local[j]>=math.ceil(BROAD_FRACTION*nh))).sum()),float(np.median(prev[mask]/nh)) if k else None,float(np.median(t[mask])) if k else None,float(np.median(b[mask])) if k else None,float(np.median(e[mask])) if k else None])
with (DEST/'threshold_summary.tsv').open('w') as f:
 w=csv.writer(f,delimiter='\t');w.writerow(['E_strictly_greater_than','fold_enrichment','route','candidates','zero_background_candidates','present_in_at_least_90pct_haps','locally_pass_in_at_least_90pct_haps','median_prevalence_fraction','median_target_count','median_background_count','median_pooled_E']);w.writerows(rows)
# Export directly from arrays: no dependence on the original selection TSV.
def seq(code):return ''.join('ACGT'[(int(code)>>shift)&3] for shift in range(30,-1,-2))
for threshold in THRESHOLDS:
 j=THRESHOLDS.index(threshold);p=(t>=POOLED_MIN)&(e>threshold);h=local[j]>0
 with (DEST/f'candidates_E_gt_{threshold}.tsv').open('w') as f:
  w=csv.writer(f,delimiter='\t');w.writerow(['kmer_id','canonical_kmer','reverse_complement','pooled_target_count','pooled_background_count','pooled_log2_enrichment','target_prevalence_haps','qualifying_haps_at_threshold','pooled_pass','hap_pass'])
  for i in np.flatnonzero(p|h):
   s=seq(codes[i]);w.writerow([f'k16_{codes[i]:08x}',s,s.translate(str.maketrans('ACGT','TGCA'))[::-1],t[i],b[i],e[i],prev[i],local[j,i],int(p[i]),int(h[i])])
with (DEST/f'candidates_E_gt_{SELECTED}.tsv').open() as f, (DEST/f'candidates_E_gt_{SELECTED}_pooled.tsv').open('w') as a, (DEST/f'candidates_E_gt_{SELECTED}_pooled_broad.tsv').open('w') as bfile:
 reader=csv.DictReader(f,delimiter='\t')
 wa=csv.DictWriter(a,fieldnames=reader.fieldnames,delimiter='\t');wb=csv.DictWriter(bfile,fieldnames=reader.fieldnames,delimiter='\t')
 wa.writeheader();wb.writeheader()
 for row in reader:
  if row['pooled_pass']=='1':
   wa.writerow(row)
   if int(row['qualifying_haps_at_threshold'])>=math.ceil(BROAD_FRACTION*nh):wb.writerow(row)
if BROAD_FRACTION==.9:
 import shutil
 shutil.copyfile(DEST/f'candidates_E_gt_{SELECTED}_pooled_broad.tsv',DEST/f'candidates_E_gt_{SELECTED}_pooled_broad90.tsv')
assert all(rows[i][3]>=rows[i+3][3] for i in range(len(rows)-3) if rows[i][2]!='hap_only')
(DEST/'validation.json').write_text(json.dumps({'status':'PASS','threshold_operator':'>','thresholds':THRESHOLDS,'selected_E':SELECTED,'nonreference_haps':nh,'preliminary_candidates':n,'epsilon':1e-9,'pooled_min_count':POOLED_MIN,'hap_min_count':HAP_MIN,'broad_hap_fraction':BROAD_FRACTION,'monotonicity_checked':True},indent=2)+'\n')
print('SWEEP_COMPLETE',flush=True)
