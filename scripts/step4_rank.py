#!/usr/bin/env python3
"""Integrate exact-match evidence, quality proxies and candidate categories."""
import csv,gzip,json,math,hashlib
from collections import defaultdict
from pathlib import Path
import numpy as np
from pipeline_schema import STEP2_METRIC_COLUMNS,require_columns
ROOT=Path(__file__).resolve().parents[1];S1=ROOT/'step1';S2=ROOT/'step2';S3=ROOT/'step3';OUT=ROOT/'step4'
for d in ['ranked','shortlists','mismatch','plots','logs']:(OUT/d).mkdir(parents=True,exist_ok=True)
def read(path):
 with (gzip.open(path,'rt') if str(path).endswith('.gz') else open(path)) as f:return list(csv.DictReader(f,delimiter='\t'))
def write(path,rows,fields=None):
 if fields is None:fields=list(rows[0]) if rows else []
 with (gzip.open(path,'wt') if str(path).endswith('.gz') else open(path,'w')) as f:
  w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',extrasaction='ignore');w.writeheader();w.writerows(rows)
def rev(s):return s.translate(str.maketrans('ACGT','TGCA'))[::-1]
def q(seq):
 gc=(seq.count('G')+seq.count('C'))/16
 # Direct run computation and stem proxy: complementary arms separated by >=3 nt.
 runs=max(len(list(g)) for _,g in __import__('itertools').groupby(seq));stem=0
 for left in range(16):
  for right in range(left+4,16):
   n=0
   while left+n<right-n and right-n<16 and seq[left+n].translate(str.maketrans('ACGT','TGCA'))==seq[right-n]:n+=1
   stem=max(stem,n)
 freq=[seq.count(b)/16 for b in 'ACGT'];entropy=-sum(x*math.log2(x) for x in freq if x)
 gcscore=max(0,1-abs(gc-.5)/.25);runscore=max(0,1-(runs-3)/4);stemscore=max(0,1-(stem-2)/4);return gc*100,runs,entropy,stem,gcscore*.35+runscore*.2+min(1,entropy/1.8)*.25+stemscore*.2
def pct(v):
 order=np.argsort(v,kind='stable');r=np.empty(len(v));r[order]=np.arange(len(v));return r/max(1,len(v)-1)
def main():
 base=read(S2/'candidate_metrics.tsv.gz');require_columns(base[0].keys() if base else [],STEP2_METRIC_COLUMNS,S2/'candidate_metrics.tsv.gz');idx={r['kmer_id']:i for i,r in enumerate(base)};codes=np.fromfile(S2/'candidate_codes.u32',dtype='<u4')
 with np.load(S1/'pooled_counts_and_scores.npz') as z:
  i=np.searchsorted(z['codes'],codes);assert np.array_equal(z['codes'][i],codes)
  bg=z['counts'][i,3,:].sum(axis=1);unionE=z['log2_enrichment'][i]
  # Use the denominator saved by STEP1.  A literal here silently made rankings
  # depend on an older STEP0 run.
  background_valid_starts=int(z['valid_starts'][3])
 categories=defaultdict(list);details={}
 for family in ['HSat2','HSat3']:
  for r in read(S3/'candidates'/f'{family}_pan_exact.tsv'):categories[r['kmer_id']].append(f'pan_{family}');details[(r['kmer_id'],f'pan_{family}')]=r
  for p in sorted((S3/'candidates').glob(f'{family}_chr*_G*_markers.tsv')):
   for r in read(p):categories[r['kmer_id']].append(f"regional_{family}_{r['chromosome']}_{r['group']}_{r['category']}");details[(r['kmer_id'],categories[r['kmer_id']][-1])]=r
  for r in read(S3/'candidates'/f'{family}_chromosome_exact.tsv'):categories[r['kmer_id']].append(f"chromosome_{family}_{r['dominant_chromosome']}");details[(r['kmer_id'],categories[r['kmer_id']][-1])]=r
 selected=[]
 for kid,cats in categories.items():
  r=base[idx[kid]].copy();seq=r['canonical_kmer'];gc,run,entropy,stem,quality=q(seq);h2=int(r['HSat2_count']);h3=int(r['HSat3_count']);family='HSat2' if h2>=h3 else 'HSat3';fi=0 if family=='HSat2' else 1
  target=max(h2,h3);other=min(h2,h3);familyE=abs(float(r['log2_HSat2_vs_HSat3_density']));familyfrac=max(float(r['HSat2_count_fraction']),1-float(r['HSat2_count_fraction']));backgroundE=math.log2((float(r[family+'_density_per_M'])/1e6+1e-9)/(bg[idx[kid]]/background_valid_starts+1e-9))
  selected.append({'kmer_id':kid,'canonical_kmer':seq,'reverse_complement':rev(seq),'target_family':family,'categories':';'.join(sorted(cats)),'target_count':target,'other_family_count':other,'family_count_fraction':familyfrac,'family_log2_enrichment':familyE,'background_count':int(bg[idx[kid]]),'family_background_log2_enrichment':backgroundE,'union_background_log2_enrichment':float(unionE[idx[kid]]),'target_prevalence_haps':int(r[family+'_prevalence_haps']),'union_prevalence_haps':int(r['union_prevalence_haps']),'GC_percent':gc,'max_homopolymer':run,'base_entropy_bits':entropy,'self_complementary_stem_proxy':stem,'sequence_quality_score':quality})
 # Scores are empirical percentiles within the integrated exact candidate universe.
 for family in ['HSat2','HSat3']:
  rows=[r for r in selected if r['target_family']==family];ab=np.log1p([r['target_count'] for r in rows]);fs=np.array([r['family_log2_enrichment'] for r in rows]);bs=np.array([r['family_background_log2_enrichment'] for r in rows]);hp=np.array([r['target_prevalence_haps'] for r in rows]);bga=-np.log1p([r['background_count'] for r in rows]);qq=np.array([r['sequence_quality_score'] for r in rows])
  components=np.vstack([pct(np.asarray(ab)),pct(fs),pct(bs),pct(hp),pct(bga),pct(qq)]).T
  for r,c in zip(rows,components):r.update(target_score=c[0],family_score=c[1],background_enrichment_score=c[2],hap_coverage_score=c[3],background_absolute_score=c[4],quality_score=c[5],exact_composite_score=float(c@np.array([.20,.20,.20,.15,.10,.15])))
 selected.sort(key=lambda r:(-r['exact_composite_score'],r['kmer_id']));write(OUT/'ranked/all_exact_ranked.tsv.gz',selected)
 # Within-category rank plus diversity selection. Similar means <=2 substitutions or shared 12-mer (including RC).
 shortlist=[]
 def similar(a,b):return min(sum(x!=y for x,y in zip(a,b)),sum(x!=y for x,y in zip(a,rev(b))))<=2 or any(a[i:i+12] in b or a[i:i+12] in rev(b) for i in range(5))
 for family in ['HSat2','HSat3']:
  catnames=sorted({c for r in selected if r['target_family']==family for c in r['categories'].split(';')})
  for cat in catnames:
   pool=[r for r in selected if cat in r['categories'].split(';')];pool.sort(key=lambda r:(-r['exact_composite_score'],r['kmer_id']));chosen=[]
   cap=20 if cat.startswith('pan_') else 5
   for r in pool:
    if r['GC_percent']<25 or r['GC_percent']>75 or r['max_homopolymer']>5 or r['self_complementary_stem_proxy']>5:continue
    if any(similar(r['canonical_kmer'],x['canonical_kmer']) for x in chosen):continue
    z=r.copy();z['shortlist_category']=cat;z['category_rank']=len(chosen)+1;chosen.append(z)
    if len(chosen)>=cap:break
   shortlist.extend(chosen)
 write(OUT/'shortlists/exact_diverse_shortlist.tsv',shortlist)
 # Mismatch expansion: all existing 24 plus best novel entries, balanced across category types and families, cap 72.
 old=read(S3/'candidates/mismatch_queries.tsv');oldids={r['kmer_id'] for r in old};queries=[]
 for r in old:
  x=next(v for v in selected if v['kmer_id']==r['kmer_id']);queries.append({**x,'query_origin':'STEP3','query_category':r['selection_category']})
 buckets=[]
 for family in ['HSat2','HSat3']:
  for typ in ['pan','regional','chromosome']:
   buckets.append([r for r in shortlist if r['target_family']==family and r['shortlist_category'].startswith(typ) and r['kmer_id'] not in oldids])
 seen=set(oldids);perbucket=defaultdict(list)
 while len(queries)<72 and any(buckets):
  progressed=False
  for bucket in buckets:
   while bucket:
    r=bucket.pop(0);key=(r['target_family'],r['shortlist_category'].split('_')[0]);prior=perbucket[key]
    if r['kmer_id'] in seen or any(similar(r['canonical_kmer'],x) for x in prior):continue
    queries.append({**r,'query_origin':'STEP4','query_category':r['shortlist_category']});seen.add(r['kmer_id']);prior.append(r['canonical_kmer']);progressed=True;break
   if len(queries)>=72:break
  if not progressed:break
 queries.sort(key=lambda r:int(r['kmer_id'][4:],16));
 for i,r in enumerate(queries):r['query_index']=i
 fields=['query_index','kmer_id','canonical_kmer','reverse_complement','target_family','query_category','query_origin','exact_composite_score','family_background_log2_enrichment','family_log2_enrichment','family_count_fraction','target_prevalence_haps','background_count','GC_percent','max_homopolymer','self_complementary_stem_proxy']
 write(OUT/'mismatch/queries.tsv',queries,fields)
 summary={'integrated_exact_candidates':len(selected),'exact_candidates_by_family':{f:sum(r['target_family']==f for r in selected) for f in ['HSat2','HSat3']},'diverse_shortlist_entries':len(shortlist),'mismatch_query_cap':72,'mismatch_queries':len(queries),'mismatch_existing_STEP3':sum(r['query_origin']=='STEP3' for r in queries),'mismatch_new_STEP4':sum(r['query_origin']=='STEP4' for r in queries),'background_valid_starts':background_valid_starts,'score_weights':{'target_abundance':.20,'family_specificity':.20,'background_enrichment':.20,'hap_coverage':.15,'absolute_background':.10,'sequence_quality_proxy':.15},'score_normalization':'empirical percentile separately within HSat2 and HSat3 integrated universes','quality_scope':'sequence-only proxy; no chemistry, Tm, RNA structure, transcriptome or expression','validation':'PASS'}
 (OUT/'ranking_summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2),flush=True)
if __name__=='__main__':main()
