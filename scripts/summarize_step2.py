#!/usr/bin/env python3
"""Build A/B matrices and specificity summaries from chromosome-resolved shards."""
import csv,gzip,json,time
import numpy as np
from step2_pipeline import OUT,S1,MANIFEST,CHROMS,REGIONS,write_tsv,sha
from workflow_settings import settings
from pipeline_schema import STEP2_CANDIDATE_COLUMNS,STEP2_METRIC_COLUMNS,require_columns

def ratio(a,b):
 return np.divide(a,b,out=np.full(np.broadcast_shapes(np.shape(a),np.shape(b)),np.nan),where=np.asarray(b)>0)
def summarize():
 codes=np.fromfile(OUT/'candidate_codes.u32',dtype='<u4');n=len(codes);nh=len(MANIFEST);nc=len(CHROMS);ref=np.array([r['sample']=='CHM13' for r in MANIFEST]);nonref=~ref
 M=OUT/'matrices';P=OUT/'plots'
 def mmap(name,shape,dtype):return np.lib.format.open_memmap(M/name,mode='w+',dtype=dtype,shape=shape)
 A=mmap('matrix_A_counts.npy',(n,nh),'<u8');AD=mmap('matrix_A_density_per_M.npy',(n,nh),'<f8')
 B=mmap('matrix_B_counts.npy',(n,nh*2),'<u8');BD=mmap('matrix_B_density_per_M.npy',(n,nh*2),'<f8')
 availability=np.zeros((nh,nc),dtype=bool);den=np.zeros((nh,nc,4),dtype=np.uint64)
 pooled=np.zeros((n,nc,4),dtype=np.uint64);dominant=np.full((n,nh,3),-1,dtype=np.int8)
 qc=[];columns=[]
 for hi,r in enumerate(MANIFEST):
  h=r['hap'];meta=json.loads((OUT/'counts'/(h+'.json')).read_text());assert meta['step1_both_orientations_match'] and meta['step0_per_chromosome_windows_match']
  with np.load(OUT/'counts'/(h+'.npz')) as z:c=z['counts'].sum(axis=3);chrs=z['chromosomes'].tolist();d=z['valid_starts']
  ci=np.array([CHROMS.index(ch) for ch in chrs]);availability[hi,ci]=True;den[hi,ci]=d
  totals=c.sum(axis=1);A[:,hi]=totals[:,2];AD[:,hi]=ratio(totals[:,2],d[:,2].sum())*1e6
  B[:,hi*2:hi*2+2]=totals[:,:2];BD[:,hi*2:hi*2+2]=ratio(totals[:,:2],d[:,:2].sum(axis=0))*1e6
  for ri in range(3):
   arr=c[:,:,ri];mx=arr.max(axis=1);best=ci[arr.argmax(axis=1)].astype(np.int8);best[(arr==mx[:,None]).sum(axis=1)>1]=-2;best[mx==0]=-1;dominant[:,hi,ri]=best
  if not ref[hi]:pooled[:,ci,:]+=c
  residual=totals[:,2]-totals[:,0]-totals[:,1]-totals[:,3]
  qc.append([h,int(ref[hi]),len(chrs),int(totals[:,2].sum()),int(residual.sum()),int(d[:,2].sum()-d[:,0].sum()-d[:,1].sum()-d[:,3].sum())])
  for fi,family in enumerate(REGIONS[:2]):
   for ch in CHROMS:
    j=chrs.index(ch) if ch in chrs else -1
    columns.append([len(columns),h,family,ch,int(j>=0),j,fi,int(d[j,fi]) if j>=0 else 'NA',f'counts/{h}.npz'])
  if (hi+1)%25==0:print('aggregate',hi+1,nh,flush=True)
 for x in [A,AD,B,BD]:x.flush()
 np.savez_compressed(M/'availability_and_denominators.npz',haplotypes=np.array([r['hap'] for r in MANIFEST]),chromosomes=np.array(CHROMS),regions=np.array(REGIONS),available=availability,valid_starts=den,is_reference=ref)
 write_tsv(M/'matrix_B_columns.tsv',['column_index','hap','family'],((hi*2+fi,r['hap'],REGIONS[fi]) for hi,r in enumerate(MANIFEST) for fi in range(2)))
 write_tsv(M/'matrix_C_columns.tsv',['column_index','hap','family','chromosome','available','shard_chromosome_index','shard_region_index','valid_starts','shard_path'],columns)
 write_tsv(OUT/'hap_qc.tsv',['hap','is_reference','available_chromosomes','selected_union_counts','selected_junction_counts','all_junction_valid_starts'],qc)
 pooled_den=den[nonref].sum(axis=0)
 density=ratio(pooled,pooled_den[None,:,:])*1e6
 np.savez_compressed(M/'pooled_chromosome_counts.npz',codes=codes,counts=pooled,valid_starts=pooled_den,density_per_M=density,chromosomes=np.array(CHROMS),regions=np.array(REGIONS))
 totals=pooled.sum(axis=1);dt=pooled_den.sum(axis=0);fd=ratio(totals[:,:2],dt[:2]);family_share=ratio(totals[:,0],totals[:,0]+totals[:,1]);fe=np.log2((fd[:,0]+1e-9)/(fd[:,1]+1e-9))
 at=np.asarray(A[:,nonref]);adt=np.asarray(AD[:,nonref]);prev=(at>0).sum(axis=1);prev10=(at>=10).sum(axis=1);prev100=(at>=100).sum(axis=1)
 tq=np.quantile(at,[.1,.5,.9],axis=1);dq=np.nanquantile(adt,[.1,.5,.9],axis=1);mean=adt.mean(axis=1);cv=ratio(adt.std(axis=1),mean)
 localb=np.asarray(B).reshape(n,nh,2)[:,nonref,:];family_prev=(localb>0).sum(axis=1)
 localden=den[nonref,:,:2].sum(axis=1);localfd=ratio(localb,localden[None,:,:]);localfe=np.log2((localfd[:,:,0]+1e-9)/(localfd[:,:,1]+1e-9))
 family_informative=((localb.sum(axis=2)>0)&np.all(localden[None,:,:]>0,axis=2)).sum(axis=1)
 h2fav=((localfe>0)&(localb.sum(axis=2)>0)).sum(axis=1);h3fav=((localfe<0)&(localb.sum(axis=2)>0)).sum(axis=1)
 with (OUT/'candidates.tsv').open() as candidate_stream:
  candidate_reader=csv.DictReader(candidate_stream,delimiter='\t');require_columns(candidate_reader.fieldnames,STEP2_CANDIDATE_COLUMNS,OUT/'candidates.tsv');candidate_rows=list(candidate_reader)
 with gzip.open(OUT/'candidate_metrics.tsv.gz','wt') as f:
  w=csv.writer(f,delimiter='\t');w.writerow(STEP2_METRIC_COLUMNS)
  for i,row in enumerate(candidate_rows):w.writerow([row['kmer_id'],row['canonical_kmer'],row['pooled_selected_threshold'],row['pooled_and_broad_threshold'],*totals[i,:3],family_share[i],*(fd[i]*1e6),fe[i],*family_prev[i],prev[i],prev[i]/nonref.sum(),prev10[i],prev100[i],*tq[:,i],*dq[:,i],cv[i],family_informative[i],h2fav[i],h3fav[i]])
 # Chromosome shares and cross-haplotype recurrence, separately for each family/union.
 chr_results={}
 with gzip.open(OUT/'chromosome_specificity.tsv.gz','wt') as f:
  w=csv.writer(f,delimiter='\t');w.writerow(['kmer_id','region','dominant_count_chromosome','dominant_count_tied','max_count_fraction','dominant_density_chromosome','dominant_density_tied','max_density_fraction','count_entropy_bits','informative_haps_for_count_chromosome','same_unique_dominant_count_chromosome_haps','dominant_count_recurrence_fraction'])
  for ri,region in enumerate(REGIONS[:3]):
   c=pooled[:,:,ri];de=density[:,:,ri];total=c.sum(axis=1);mx=c.max(axis=1);best=c.argmax(axis=1);ties=(c==mx[:,None]).sum(axis=1)>1
   dm=np.nanmax(de,axis=1);db=np.nanargmax(de,axis=1);dti=(de==dm[:,None]).sum(axis=1)>1
   fractions=ratio(c,total[:,None]);entropy=-np.nansum(np.where(fractions>0,fractions*np.log2(np.where(fractions>0,fractions,1)),0),axis=1)
   if ri<2:present=localb[:,:,ri]>0
   else:present=at>0
   region_available=den[nonref,:,ri]>0
   informative=present&region_available[:,best].T&(region_available.sum(axis=1)[None,:]>=2)
   informative[ties|(mx==0)]=False
   support=((dominant[:,nonref,ri]==best[:,None])&informative).sum(axis=1);eligible=informative.sum(axis=1);recurrence=ratio(support,eligible)
   share=ratio(mx,total);dshare=ratio(dm,np.nansum(de,axis=1));chr_results[region]=(share,recurrence)
   for i in range(n):
    w.writerow([candidate_rows[i]['kmer_id'],region,CHROMS[best[i]] if total[i] else 'NA',int(ties[i]) if total[i] else 'NA',share[i],CHROMS[db[i]] if dm[i]>0 else 'NA',int(dti[i]) if dm[i]>0 else 'NA',dshare[i],entropy[i] if total[i] else 'NA',eligible[i],support[i],recurrence[i]])
 # Small reproducible plot tables; no clustering is performed here.
 rng=np.random.default_rng(42);sample=np.sort(rng.choice(n,min(n,40000),replace=False))
 write_tsv(P/'family_scatter.tsv',['log10_HSat2_density_plus_epsilon','log10_HSat3_density_plus_epsilon','pooled_pass'],((np.log10(fd[i,0]+1e-9),np.log10(fd[i,1]+1e-9),candidate_rows[i]['pooled_selected_threshold']) for i in sample))
 def histogram(name,values,edges):
  h,e=np.histogram(values[np.isfinite(values)],bins=edges);write_tsv(P/name,['midpoint','count'],zip((e[:-1]+e[1:])/2,h))
 histogram('family_share_hist.tsv',family_share,np.linspace(0,1,21));histogram('prevalence_hist.tsv',prev/nonref.sum(),np.linspace(0,1,21));histogram('union_chr_share_hist.tsv',chr_results['HSat23'][0],np.linspace(0,1,21))
 write_tsv(P/'specificity_prevalence.tsv',['prevalence_fraction','max_chromosome_count_fraction','dominant_chromosome_recurrence'],((prev[i]/nonref.sum(),chr_results['HSat23'][0][i],chr_results['HSat23'][1][i]) for i in sample))
 pooled_mask=np.array([r['pooled_selected_threshold']=='1' for r in candidate_rows]);ids=[]
 for ri in [0,1]:
  rank=np.argsort(-fd[:,ri],kind='stable');ids.extend([int(i) for i in rank if pooled_mask[i]][:15])
 ids=list(dict.fromkeys(ids));heat=ratio(density[ids,:,2],np.nanmax(density[ids,:,2],axis=1)[:,None])
 write_tsv(P/'chromosome_heatmap.tsv',['kmer_id']+CHROMS,([candidate_rows[i]['kmer_id']]+heat[j].tolist() for j,i in enumerate(ids)))
 broad_min=int(np.ceil(settings()['selection']['broad_hap_fraction']*int(nonref.sum())))
 summary={'candidate_count':n,'sequence_sets':nh,'nonreference_haps':int(nonref.sum()),'broad_haps_min':broad_min,'reference':'chm13v2.0','observed_hap_chromosome_pairs':int(availability.sum()),'pooled_selected_candidates':int(pooled_mask.sum()),'pooled_and_broad_threshold_candidates':sum(int(r['pooled_and_broad_threshold']) for r in candidate_rows),'ambiguous_valid_starts':int(den[:,:,3].sum()),'HSat2_count_fraction_ge_0.9':int((family_share>=.9).sum()),'HSat3_count_fraction_ge_0.9':int((family_share<=.1).sum()),'union_max_chr_count_fraction_ge_0.9':int((chr_results['HSat23'][0]>=.9).sum()),'union_prevalence_ge_broad_threshold':int((prev>=broad_min).sum()),'matrix_A_shape':[n,nh],'matrix_B_shape':[n,nh*2],'matrix_C_family_columns_including_missing':len(columns),'step1_per_hap_agreement':True,'counts_dtype':'uint64','density_unit':'counts per million valid 16-mer starts'}
 (OUT/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2),flush=True)
if __name__=='__main__':summarize()
