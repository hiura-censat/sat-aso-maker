#!/usr/bin/env python3
"""HSat2-selective exact-match shortlist, not finalized ASO chemistry/strand."""
import csv,gzip,json,hashlib
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];S1=ROOT/'step1';S2=ROOT/'step2';OUT=S2/'hsat2_aso_candidates';OUT.mkdir(exist_ok=True)
with gzip.open(S2/'candidate_metrics.tsv.gz','rt') as f:
 rows=[r for r in csv.DictReader(f,delimiter='\t') if r['pooled_and_broad90_E_gt_10']=='1' and float(r['HSat2_count_fraction'])>=.999 and float(r['log2_HSat2_vs_HSat3_density'])>10 and int(r['HSat2_prevalence_haps'])>=516]
codes=np.array([int(r['kmer_id'][4:],16) for r in rows],dtype=np.uint32)
allcodes=np.fromfile(S1/'preliminary_candidates.u32',dtype='<u4');i1=np.searchsorted(allcodes,codes);assert np.array_equal(allcodes[i1],codes)
c2=np.fromfile(S2/'candidate_codes.u32',dtype='<u4');i2=np.searchsorted(c2,codes);assert np.array_equal(c2[i2],codes)
with np.load(S1/'pooled_counts_and_scores.npz') as z:bg=z['counts'][i1,3,:].sum(axis=1);bd=z['background_density'][i1]
B=np.load(S2/'matrices/matrix_B_counts.npy',mmap_mode='r')[i2].reshape(len(rows),574,2)
BD=np.load(S2/'matrices/matrix_B_density_per_M.npy',mmap_mode='r')[i2].reshape(len(rows),574,2)/1e6
haps=list(csv.DictReader((S2/'haplotypes.tsv').open(),delimiter='\t'));nonref=np.array([r['is_reference']=='0' for r in haps]);assert nonref.sum()==573
local_bg=np.zeros(len(rows),dtype=np.uint16);local_both=local_bg.copy()
for hi,h in enumerate(haps):
 if h['is_reference']=='1':continue
 with np.load(S1/'enrichment_by_hap'/(h['hap']+'.npz')) as z:bgd=z['background_density'][i1]
 ehb=np.log2((BD[:,hi,0]+1e-9)/(bgd+1e-9));eh3=np.log2((BD[:,hi,0]+1e-9)/(BD[:,hi,1]+1e-9))
 denom=B[:,hi,:].sum(axis=1);frac=np.divide(B[:,hi,0],denom,out=np.zeros(len(rows)),where=denom>0)
 passed=(B[:,hi,0]>=10)&(ehb>10);local_bg+=passed.astype(np.uint16);local_both+=(passed&(eh3>10)&(frac>=.999)).astype(np.uint16)
 if (hi+1)%100==0:print('hap',hi+1,574,flush=True)
quant=np.quantile(B[:,nonref,0],[.1,.5,.9],axis=1)
def reverse(s):return s.translate(str.maketrans('ACGT','TGCA'))[::-1]
for i,r in enumerate(rows):
 seq=r['canonical_kmer'];r.update(reverse_complement=reverse(seq),pooled_background_count=int(bg[i]),HSat2_vs_background_E=float(np.log2((float(r['HSat2_density_per_M'])/1e6+1e-9)/(bd[i]+1e-9))),HSat2_E_gt10_count_ge10_haps=int(local_bg[i]),HSat2_enrichment_and_strict_family_haps=int(local_both[i]),HSat2_count_p10=quant[0,i],HSat2_count_median=quant[1,i],HSat2_count_p90=quant[2,i],GC_percent=100*(seq.count('G')+seq.count('C'))/16)
rows=[r for r in rows if r['HSat2_vs_background_E']>10]
rows.sort(key=lambda r:(-r['HSat2_vs_background_E'],r['kmer_id']))
broad=[r for r in rows if r['HSat2_E_gt10_count_ge10_haps']>=516]
strict=[r for r in broad if r['HSat2_enrichment_and_strict_family_haps']>=516]
for name,data in [('HSat2_selective_pooled.tsv',rows),('HSat2_selective_broad90.tsv',broad),('HSat2_selective_strict_broad90.tsv',strict)]:
 with (OUT/name).open('w') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(data)
# Independent identity and ordering checks, plus monotone subset membership.
assert set(r['kmer_id'] for r in strict)<=set(r['kmer_id'] for r in broad)<=set(r['kmer_id'] for r in rows)
for i,r in enumerate(rows):assert len(r['canonical_kmer'])==16 and r['canonical_kmer']<=r['reverse_complement'] and 0<=r['HSat2_enrichment_and_strict_family_haps']<=r['HSat2_E_gt10_count_ge10_haps']<=573
summary={'pooled_selective':len(rows),'broad_HSat2_enrichment':len(broad),'broad_HSat2_enrichment_and_family_specificity':len(strict),'criteria':{'HSat2_count_fraction_min':.999,'HSat2_vs_HSat3_E_strictly_greater_than':10,'HSat2_vs_background_E_strictly_greater_than':10,'local_HSat2_count_min':10,'broad_haps_min':516,'total_nonreference_haps':573,'epsilon':1e-9,'initial_universe':'STEP2 pooled_and_broad90_E_gt_10 candidates'},'rank':'descending pooled HSAT2 vs background E; no diversity or chemistry filter','validation':'PASS','source_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [S2/'candidate_metrics.tsv.gz',S2/'COMPLETE.json',S1/'COMPLETE.json']}}
(OUT/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2),flush=True)
for r in strict[:10]:print(r['canonical_kmer'],r['reverse_complement'],r['HSat2_count'],r['HSat3_count'],r['pooled_background_count'],round(r['HSat2_vs_background_E'],3),r['HSat2_E_gt10_count_ge10_haps'],r['HSat2_enrichment_and_strict_family_haps'],flush=True)
