#!/usr/bin/env python3
import csv,json
import numpy as np
from step3_prepare import ROOT,OUT,tsv
from step2_pipeline import MANIFEST,STATS,CHROMS
from workflow_settings import settings
M=OUT/'mismatch'
def summarize():
 queries=list(csv.DictReader((OUT/'candidates/mismatch_queries.tsv').open(),delimiter='\t'));nq=len(queries);nh=len(MANIFEST)
 allcounts=np.zeros((nq,nh,4,3),dtype=np.uint64);den=np.zeros((nh,4),dtype=np.uint64)
 positions=np.zeros((nq,nh,4,3,16));reference=np.array([r['sample']=='CHM13' for r in MANIFEST]);checked=[]
 broad_min=int(np.ceil(settings()['selection']['broad_hap_fraction']*int((~reference).sum())))
 chromosome_counts=np.zeros((nq,nh,len(CHROMS),4,3),dtype=np.uint64);chromosome_den=np.zeros((nh,len(CHROMS),4),dtype=np.uint64);available=np.zeros((nh,len(CHROMS)),dtype=bool)
 for hi,r in enumerate(MANIFEST):
  h=r['hap'];qc=json.loads((M/'counts'/(h+'.json')).read_text());assert qc['zero_mismatch_matches_step1'] and qc['all_valid_windows_match_step0'] and qc['mismatch_position_accounting']
  with np.load(M/'counts'/(h+'.npz')) as z:c=z['counts'];d=z['valid_starts'];ch=z['chromosomes'];p=z['position_counts']
  assert c.shape==(nq,len(ch),4,3) and d.shape==(len(ch),4);assert np.all(c.sum(axis=3)<=d[None,:,:])
  for distance in range(3):np.testing.assert_array_equal(p[:,:,distance,:].sum(axis=2),c[:,:,:,distance].sum(axis=1)*distance)
  allcounts[:,hi]=c.sum(axis=1);den[hi]=d.sum(axis=0);positions[:,hi]=p
  ci=[CHROMS.index(str(s)) for s in ch];chromosome_counts[:,hi,ci,:,:]=c;chromosome_den[hi,ci]=d;available[hi,ci]=True
  checked.append(h)
 pool=allcounts[:,~reference].sum(axis=1);pd=den[~reference].sum(axis=0);cum=pool.cumsum(axis=2);lc=allcounts.cumsum(axis=3)
 rows=[];decisions=[]
 for qi,q in enumerate(queries):
  fi=['HSat2','HSat3'].index(q['family']);other=1-fi;flags=[]
  for radius in range(3):
   t=cum[qi,fi,radius];o=cum[qi,other,radius];b=cum[qi,3,radius];e=np.log2((t/pd[fi]+1e-9)/(b/pd[3]+1e-9));ef=np.log2((t/pd[fi]+1e-9)/(o/pd[other]+1e-9));fraction=t/(t+o) if t+o else np.nan
   with np.errstate(divide='ignore',invalid='ignore'):
    le=np.log2((lc[qi,:,fi,radius]/den[:,fi]+1e-9)/(lc[qi,:,3,radius]/den[:,3]+1e-9));lf=np.log2((lc[qi,:,fi,radius]/den[:,fi]+1e-9)/(lc[qi,:,other,radius]/den[:,other]+1e-9));lfrac=lc[qi,:,fi,radius]/(lc[qi,:,fi,radius]+lc[qi,:,other,radius])
   local=(lc[qi,:,fi,radius]>=10)&(le>10)&(lf>10)&(lfrac>=.999)&~reference
   broad=int(local.sum());passed=bool(e>10 and ef>10 and fraction>=.999 and broad>=broad_min);flags.append(passed)
   rows.append([qi,q['kmer_id'],q['canonical_kmer'],q['family'],q['selection_category'],radius,int(t),int(o),int(b),int(pool[qi,fi,radius]),int(pool[qi,other,radius]),int(pool[qi,3,radius]),float(e),float(ef),float(fraction),broad,int(passed)])
  decisions.append([q['kmer_id'],q['family'],q['selection_category'],*map(int,flags)])
 tsv(M/'pooled_mismatch_scores.tsv',['query_index','kmer_id','canonical_kmer','family','selection_category','max_mismatches','target_count_le_radius','other_family_count_le_radius','background_count_le_radius','target_count_exact_distance','other_family_count_exact_distance','background_count_exact_distance','target_vs_background_E','target_vs_other_family_E','target_family_count_fraction','strict_passing_nonreference_haps','passes_strict_pan_criteria'],rows)
 tsv(M/'candidate_robustness.tsv',['kmer_id','family','selection_category','strict_pan_pass_0mm','strict_pan_pass_le1mm','strict_pan_pass_le2mm'],decisions)
 # Chromosome-group contrasts compare only haplotypes measurable on that chromosome.
 group_rows=[];group_robust=[];case_summaries=json.loads((OUT/'clustering_summary.json').read_text())
 for case in case_summaries:
  if not case['supported_partition']:continue
  family=case['family'];chrom=case['chromosome'];fi=['HSat2','HSat3'].index(family);ci=CHROMS.index(chrom)
  with np.load(OUT/'clustering'/case['case']/'model.npz') as z:hi=z['hap_indices'];labels=z['labels']
  assert np.all(chromosome_den[hi,ci,fi]>0) and not reference[hi].any()
  for qi,q in enumerate(queries):
   if q['family']!=family:continue
   ct=chromosome_counts[qi,hi,ci,fi,:].cumsum(axis=1);density=ct/chromosome_den[hi,ci,fi,None]*1e6
   for g in np.unique(labels):
    inside=labels==g
    for radius in range(3):
     pin=float((ct[inside,radius]>=10).mean());pout=float((ct[~inside,radius]>=10).mean())
     a=float(np.median(density[inside,radius]));b=float(np.median(density[~inside,radius]))
     score=float(np.log2((a+.001)/(b+.001)))
     group_rows.append([qi,q['kmer_id'],family,chrom,f'G{g}',radius,pin,pout,a,b,score,int(pin>=.8 and pout<=.2 and score>=2)])
 tsv(M/'group_contrast_by_mismatch.tsv',['query_index','kmer_id','family','chromosome','group','max_mismatches','in_group_count_ge10_fraction','out_group_count_ge10_fraction','in_group_median_density_per_M','out_group_median_density_per_M','log2_median_density_ratio','passes_group_presence_contrast'],group_rows)
 index={(int(r[0]),r[3],r[4],int(r[5])):r for r in group_rows}
 score_index={(int(r[0]),int(r[5])):r for r in rows}
 for qi,q in enumerate(queries):
  category=q['selection_category']
  if not category.startswith('chr') or category.startswith('chromosome_'):continue
  chrom,group,*_=category.split('_')
  exclusive=category.endswith('_specific')
  for radius in range(3):
   r=index.get((qi,chrom,group,radius))
   if r is None:continue
   sr=score_index[(qi,radius)];family_ok=sr[12]>10 and sr[13]>6 and sr[14]>=.99
   regional=r[6]>=.8 and r[10]>=2 and (not exclusive or r[7]<=.2)
   group_robust.append([qi,q['kmer_id'],q['family'],chrom,group,category,radius,int(family_ok),int(regional),int(family_ok and regional)])
 tsv(M/'selected_group_robustness.tsv',['query_index','kmer_id','family','chromosome','group','selection_category','max_mismatches','passes_original_family_background_screen','retains_chromosome_group_contrast','passes_family_and_chromosome_group_screen'],group_robust)
 # Chromosome specificity is conditional on measured annotated regions.
 pooled_chr=chromosome_counts[:,~reference].sum(axis=1).cumsum(axis=3);local_chr=chromosome_counts.cumsum(axis=4);chrrows=[]
 for qi,q in enumerate(queries):
  fi=['HSat2','HSat3'].index(q['family']);measured=chromosome_den[:,:,fi]>0
  for radius in range(3):
   c=pooled_chr[qi,:,fi,radius];total=c.sum();best=int(c.argmax());tied=int((c==c[best]).sum())>1;lc=local_chr[qi,:,:,fi,radius];mx=lc.max(axis=1);lb=lc.argmax(axis=1);lt=(lc==mx[:,None]).sum(axis=1)>1
   eligible=~reference&(lc.sum(axis=1)>0)&measured[:,best]&(measured.sum(axis=1)>=2)
   if tied:eligible[:]=False
   support=int((eligible&~lt&(lb==best)).sum());count=int(eligible.sum());fraction=float(c[best]/total) if total else np.nan;recurrence=support/count if count else np.nan
   chrrows.append([qi,q['kmer_id'],q['family'],q['selection_category'],radius,CHROMS[best] if total else 'NA',int(tied),fraction,count,support,recurrence,int(fraction>=.9 and recurrence>=.8 and count>=100)])
 tsv(M/'chromosome_robustness.tsv',['query_index','kmer_id','family','selection_category','max_mismatches','dominant_chromosome','dominant_tied','max_target_count_fraction','informative_haps','same_unique_dominant_chromosome_haps','recurrence_fraction','passes_chromosome_concentration_screen'],chrrows)
 pp=positions[:,~reference].sum(axis=1)
 tsv(M/'mismatch_position_counts.tsv',['query_index','kmer_id','region','exact_distance','query_position_1based','weighted_count'],((qi,q['kmer_id'],region,d,pos+1,float(pp[qi,ri,d,pos])) for qi,q in enumerate(queries) for ri,region in enumerate(['HSat2','HSat3','HSat23','background']) for d in [1,2] for pos in range(16)))
 for d in range(3):np.testing.assert_array_equal(pp[:,:,d,:].sum(axis=2),pool[:,:,d]*d)
 np.savez_compressed(M/'hap_counts_and_positions.npz',counts=allcounts,valid_starts=den,position_counts=positions,chromosome_counts=chromosome_counts,chromosome_valid_starts=chromosome_den,chromosome_available=available,chromosomes=np.array(CHROMS),query_codes=np.array([int(q['kmer_id'][4:],16) for q in queries],dtype=np.uint32),haplotypes=np.array([r['hap'] for r in MANIFEST]),is_reference=reference)
 summary={'status':'PASS','queries':nq,'sequence_sets':nh,'nonreference_haps':int((~reference).sum()),'broad_haps_min':broad_min,'all_distance_zero_counts_match_step1':True,'all_valid_windows_match_step0':True,'position_count_identity_verified':True,'strict_pan_pass_counts_by_radius':[sum(int(r[3+i]) for r in decisions) for i in range(3)],'method':'exhaustive Hamming radius <=2 against both orientations, canonical window counted once per query; no indels','scope':'eligible T2T chromosomes, existing HSat2/3 annotations and their background complement'}
 summary['pan_selected_queries']=sum(q['selection_category']=='pan' for q in queries);summary['pan_selected_strict_pass_counts_by_radius']=[sum(int(r[3+i]) for r in decisions if r[2]=='pan') for i in range(3)]
 summary['selected_group_queries']=sum(q['selection_category'].startswith('chr') and not q['selection_category'].startswith('chromosome_') for q in queries);summary['group_family_and_chromosome_pass_by_radius']=[sum(r[9] for r in group_robust if r[6]==i) for i in range(3)]
 (M/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2),flush=True)
if __name__=='__main__':summarize()
