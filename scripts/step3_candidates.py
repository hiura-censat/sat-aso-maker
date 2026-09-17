#!/usr/bin/env python3
import csv,gzip,json
import numpy as np
from step3_prepare import ROOT,OUT,S2,tsv
from workflow_settings import settings

def sequence(code):return ''.join('ACGT'[(int(code)>>i)&3] for i in range(30,-1,-2))
def rc(s):return s.translate(str.maketrans('ACGT','TGCA'))[::-1]
def sufficiently_distinct(s,existing):
 for t in existing:
  for u in [t,rc(t)]:
   if sum(a!=b for a,b in zip(s,u))<=2:return False
   if any(s[i:i+12] in u for i in range(5)):return False
 return True
def candidates():
 codes=np.fromfile(S2/'candidate_codes.u32',dtype='<u4');n=len(codes)
 with np.load(OUT/'data/metadata.npz') as z:ref=z['is_reference'];den=z['valid_starts'];haps=z['haplotypes'];chromosomes=list(z['chromosomes'])
 B=np.load(S2/'matrices/matrix_B_counts.npy',mmap_mode='r')
 chrom_counts=np.load(OUT/'data/chromosome_feature_counts.npy',mmap_mode='r')
 feature_pool=np.load(OUT/'data/chromosome_feature_indices.npy')
 cluster_summaries=json.loads((OUT/'clustering_summary.json').read_text())
 with gzip.open(S2/'candidate_metrics.tsv.gz','rt') as f:metrics=list(csv.DictReader(f,delimiter='\t'))
 with np.load(ROOT/'step1/pooled_counts_and_scores.npz') as z:
  ix=np.searchsorted(z['codes'],codes);assert np.array_equal(z['codes'][ix],codes);background=z['background_density'][ix];background_count=z['counts'][ix,3,:].sum(axis=1)
 family_counts=np.array([[int(r['HSat2_count']),int(r['HSat3_count'])] for r in metrics]);fd=family_counts/den[~ref,:,:2].sum(axis=(0,1))[None,:]
 es=np.log2((fd+1e-9)/(background[:,None]+1e-9));fe=np.log2((fd+1e-9)/(fd[:,::-1]+1e-9));shares=np.divide(family_counts,family_counts.sum(axis=1)[:,None],out=np.zeros_like(fd),where=family_counts.sum(axis=1)[:,None]>0)
 nonreference_haps=int((~ref).sum());broad_fraction=settings()['selection']['broad_hap_fraction'];broad_min=int(np.ceil(broad_fraction*nonreference_haps))
 pools=[];summary=[];all_group_rows=[]
 for fi,family in enumerate(['HSat2','HSat3']):
  C=np.asarray(B[:,fi::2][:,~ref]);broad=(C>=10).sum(axis=1)>=broad_min;family_ok=(shares[:,fi]>=.99)&(fe[:,fi]>6)&(es[:,fi]>10)
  strict=(shares[:,fi]>=.999)&(fe[:,fi]>10)&(es[:,fi]>10)&broad
  ids=np.flatnonzero(strict);ids=ids[np.argsort(-es[ids,fi],kind='stable')]
  tsv(OUT/'candidates'/(family+'_pan_exact.tsv'),['kmer_id','canonical_kmer','family_background_E','family_count_fraction','other_family_E','count_ge10_hap_fraction','pooled_background_count'],((f'k16_{codes[i]:08x}',sequence(codes[i]),es[i,fi],shares[i,fi],fe[i,fi],float((C[i]>=10).mean()),int(background_count[i])) for i in ids))
  pools.append((family,'pan',ids.tolist(),fi))
  for case in cluster_summaries:
   if case['family']!=family or not case['supported_partition']:continue
   chrom=case['chromosome'];ci=chromosomes.index(chrom)
   with np.load(OUT/'clustering'/case['case']/'model.npz') as z:hi=z['hap_indices'];labels=z['labels']
   cd=den[hi,ci,fi];c=np.asarray(chrom_counts[hi,ci,fi,:]).T;density=c/cd[None,:]*1e6;index=feature_pool[ci,fi]
   for g in np.unique(labels):
    inside=labels==g;pin=(c[:,inside]>=10).mean(axis=1);pout=(c[:,~inside]>=10).mean(axis=1);medin=np.median(density[:,inside],axis=1);medout=np.median(density[:,~inside],axis=1);score=np.log2((medin+.001)/(medout+.001))
    ok=family_ok[index]&(pin>=.8)&(score>=2);selected=np.flatnonzero(ok);selected=selected[np.argsort(-score[selected],kind='stable')]
    rows=[]
    for j in selected:
     i=index[j];category='group_specific_exact' if pout[j]<=.2 else 'group_enriched_exact';rows.append([f'k16_{codes[i]:08x}',sequence(codes[i]),family,chrom,f'G{g}',category,int(inside.sum()),int((~inside).sum()),pin[j],pout[j],medin[j],medout[j],score[j],es[i,fi],shares[i,fi]])
    tsv(OUT/'candidates'/f'{family}_{chrom}_G{g}_markers.tsv',['kmer_id','canonical_kmer','family','chromosome','group','category','in_group_haps','out_group_haps','in_group_count_ge10_fraction','out_group_count_ge10_fraction','in_group_median_density_per_M','out_group_median_density_per_M','log2_median_density_ratio','family_background_E','family_count_fraction'],rows)
    specific=[index[j] for j in selected if pout[j]<=.2];relaxed=[index[j] for j in selected if pout[j]>.2];category=f'{chrom}_G{g}_'+('specific' if specific else 'enriched');pools.append((family,category,specific if specific else relaxed,fi))
    all_group_rows.append([family,chrom,f'G{g}',int(inside.sum()),len(specific),len(relaxed)])
  summary.append({'family':family,'pan_exact_candidates':len(ids),'chromosome_models':sum(r['family']==family and r['supported_partition'] for r in cluster_summaries)})
 # Retain strong chromosome-concentrated candidates for mismatch comparison.
 chrrows=[[],[]]
 with gzip.open(S2/'chromosome_specificity.tsv.gz','rt') as f:
  for r in csv.DictReader(f,delimiter='\t'):
   if r['region'] not in ['HSat2','HSat3']:continue
   fi=['HSat2','HSat3'].index(r['region']);code=int(r['kmer_id'][4:],16);i=int(np.searchsorted(codes,code))
   if float(r['max_count_fraction'])>=.9 and float(r['dominant_count_recurrence_fraction'])>=.8 and int(r['informative_haps_for_count_chromosome'])>=100 and shares[i,fi]>=.999 and fe[i,fi]>10 and es[i,fi]>10:chrrows[fi].append((i,r['dominant_count_chromosome'],float(r['max_count_fraction']),float(r['dominant_count_recurrence_fraction'])))
 for fi,family in enumerate(['HSat2','HSat3']):
  chrrows[fi].sort(key=lambda v:-es[v[0],fi]);tsv(OUT/'candidates'/(family+'_chromosome_exact.tsv'),['kmer_id','canonical_kmer','dominant_chromosome','count_fraction','recurrence','family_background_E'],((f'k16_{codes[i]:08x}',sequence(codes[i]),ch,a,b,es[i,fi]) for i,ch,a,b in chrrows[fi]))
  # Preserve the intended chromosome in the STEP3 mismatch provenance.  Each
  # chromosome is a separate category, just like the regional group pools.
  for chrom in chromosomes:
   ids=[i for i,ch,_,_ in chrrows[fi] if ch==chrom]
   if ids:pools.append((family,f'chromosome_{chrom}',ids,fi))
 # Round-robin categories, at most 24 distinct queries; first round includes each nonempty pool.
 chosen=[];seen=set();provenance={};perpool={}
 for roundno in range(3):
  for family,category,ids,fi in pools:
   if len(chosen)>=24:break
   key=(family,category);previous=perpool.setdefault(key,[])
   for i in ids:
    if i in seen:continue
    s=sequence(codes[i])
    if not sufficiently_distinct(s,previous):continue
    chosen.append(i);seen.add(i);previous.append(s);provenance[i]=(family,category);break
 chosen.sort(key=lambda i:int(codes[i]))
 tsv(OUT/'candidates/mismatch_queries.tsv',['query_index','kmer_id','canonical_kmer','reverse_complement','family','selection_category','exact_family_background_E'],((j,f'k16_{codes[i]:08x}',sequence(codes[i]),rc(sequence(codes[i])),*provenance[i],es[i,['HSat2','HSat3'].index(provenance[i][0])]) for j,i in enumerate(chosen)))
 tsv(OUT/'candidates/group_marker_summary.tsv',['family','chromosome','group','haps','group_specific_exact_candidates','group_enriched_exact_candidates'],all_group_rows)
 (OUT/'candidates/summary.json').write_text(json.dumps({'families':summary,'nonreference_haps':nonreference_haps,'broad_hap_fraction':broad_fraction,'broad_haps_min':broad_min,'mismatch_query_count':len(chosen),'query_cap':24,'note':'same-data exploratory effect sizes, not independent validation; no p-values'},indent=2)+'\n');print('CANDIDATES_COMPLETE',len(chosen),flush=True)
if __name__=='__main__':candidates()
