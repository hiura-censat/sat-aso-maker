#!/usr/bin/env python3
import csv,json,subprocess
from pathlib import Path
import numpy as np
from step3_prepare import ROOT,OUT,S2,tsv

def ari(a,b):
 _,a=np.unique(a,return_inverse=True);_,b=np.unique(b,return_inverse=True);t=np.zeros((a.max()+1,b.max()+1),dtype=np.int64);np.add.at(t,(a,b),1)
 c2=lambda x:np.sum(x*(x-1)/2);n=len(a);expect=c2(t.sum(axis=0))*c2(t.sum(axis=1))/(n*(n-1)/2);den=(c2(t.sum(axis=0))+c2(t.sum(axis=1)))/2-expect
 return float((c2(t)-expect)/den) if den else float(np.array_equal(a[:,None]==a[None,:],b[:,None]==b[None,:]))
def run_case(name,counts,den,hap_indices,feature_codes,meta,reps=0,reference=None):
 folder=OUT/'clustering'/name;folder.mkdir(exist_ok=True)
 if (folder/'model.npz').exists() and (folder/'summary.json').exists():return json.loads((folder/'summary.json').read_text())
 counts=np.asarray(counts,dtype=np.float64);x=np.log1p(counts/den[:,None]*1e6)
 prevalence=(counts>=10).mean(axis=0);variance=x.var(axis=0);ok=(prevalence>=.05)&(variance>1e-4)
 selected=np.flatnonzero(ok);selected=selected[np.argsort(-variance[selected],kind='stable')[:3000]]
 if len(selected)<2:
  summary={'case':name,'haps':len(hap_indices),'selected_k':1,'supported_partition':False,'exploratory_only':reps==0,'features':len(selected),'status':'insufficient_variable_features_in_analyzed_pool'}
  (folder/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');tsv(folder/'assignments.tsv',['hap','group','status'],((str(h),'NA',summary['status']) for h in meta['haplotypes'][hap_indices]));print('UNCLASSIFIED',name,flush=True);return summary
 pool=x[:,selected];standard=(pool-pool.mean(axis=0))/pool.std(axis=0);covered=np.zeros(len(selected),dtype=bool);representatives=[];mapping=[]
 for j in range(len(selected)):
  if covered[j]:continue
  correlations=standard.T@standard[:,j]/len(pool);members=np.flatnonzero((np.abs(correlations)>=.95)&~covered);covered[members]=True;representatives.append(j)
  mapping.extend((int(feature_codes[selected[m]]),int(feature_codes[selected[j]]),float(correlations[m])) for m in members)
  if len(representatives)>=500:break
 ix=selected[representatives];features=x[:,ix];mean=features.mean(axis=0);u,s,vh=np.linalg.svd(features-mean,full_matrices=False);scores=u*s;explained=s*s/(s*s).sum()
 haps=meta['haplotypes'][hap_indices];tsv(folder/'features.tsv',['hap']+[f'k16_{c:08x}' for c in feature_codes[ix]],([h]+features[i].tolist() for i,h in enumerate(haps)))
 tsv(folder/'feature_groups.tsv',['kmer_id','representative_kmer_id','correlation'],((f'k16_{c:08x}',f'k16_{r:08x}',v) for c,r,v in mapping))
 tsv(folder/'pca.tsv',['hap']+[f'PC{i+1}' for i in range(min(10,scores.shape[1]))],([h]+scores[i,:10].tolist() for i,h in enumerate(haps)))
 tsv(folder/'pca_explained.tsv',['PC','variance_fraction'],((i+1,v) for i,v in enumerate(explained)))
 subprocess.run(['Rscript',str(ROOT/'scripts/step3_hclust.R'),str(folder),str(reps)],check=True)
 results=list(csv.DictReader((folder/'cluster_selection.tsv').open(),delimiter='\t'));parts=list(csv.DictReader((folder/'partitions.tsv').open(),delimiter='\t'))
 eligible=[r for r in results if float(r['silhouette'])>=.25 and int(r['min_cluster_size'])>=10 and (reps==0 or (float(r['median_subsample_ARI'])>=.75 and float(r['median_worst_cluster_Jaccard'])>=.75))]
 best=max(eligible,key=lambda r:(float(r['silhouette']),-int(r['k']))) if eligible else None;k=int(best['k']) if best else 1
 labels=np.array([int(r['k'+str(k)]) for r in parts]) if k>1 else np.ones(len(haps),dtype=int)
 order=sorted(np.unique(labels),key=lambda j:scores[labels==j,0].mean());lookup={g:i+1 for i,g in enumerate(order)};labels=np.array([lookup[g] for g in labels])
 tsv(folder/'assignments.tsv',['hap','group','status'],((h,f'G{g}' if k>1 else 'unresolved','stable_partition' if reps and k>1 else 'exploratory_partition' if k>1 else 'no_supported_partition') for h,g in zip(haps,labels)))
 ref_info={}
 if reference is not None:
  rc,rd=reference;rf=np.log1p(rc[ix]/rd*1e6);rp=(rf-mean)@vh.T;centroids=np.array([features[labels==j].mean(axis=0) for j in range(1,k+1)]);dist=np.linalg.norm(centroids-rf,axis=1);nearest=int(dist.argmin()+1);train=np.linalg.norm(features[labels==nearest]-centroids[nearest-1],axis=1);ref_info={'nearest_group':f'G{nearest}' if k>1 else 'unresolved','distance':float(dist.min()),'outside_training_group_95pct_distance':bool(dist.min()>np.quantile(train,.95))}
  tsv(folder/'CHM13_projection.tsv',['reference']+[f'PC{i+1}' for i in range(min(10,len(rp)))],[['chm13v2.0']+rp[:10].tolist()])
 np.savez_compressed(folder/'model.npz',hap_indices=hap_indices,labels=labels,feature_codes=feature_codes[ix],feature_indices=ix,features=features,mean=mean,components=vh,scores=scores,explained=explained)
 summary={'case':name,'haps':len(haps),'selected_k':k,'supported_partition':bool(k>1 and reps>0),'exploratory_only':reps==0,'features':len(ix),'variable_features':int(ok.sum()),'correlation_pool_features':len(selected),'mapped_correlated_features':int(covered.sum()),'feature_subsamples':reps,'selection':best,'reference_projection':ref_info}
 (folder/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print('CLUSTER',name,len(haps),k,len(ix),flush=True);return summary
def cluster():
 with np.load(OUT/'data/metadata.npz') as z:meta={k:z[k] for k in z.files}
 C=np.load(OUT/'data/core_counts.npy',mmap_mode='r');CC=np.load(OUT/'data/chromosome_feature_counts.npy',mmap_mode='r');pool=np.load(OUT/'data/chromosome_feature_indices.npy');codes=meta['codes'];ref=meta['is_reference'];summaries=[]
 B=np.load(S2/'matrices/matrix_B_counts.npy',mmap_mode='r')
 for fi,family in enumerate(['HSat2','HSat3']):
  use=np.flatnonzero(~ref&meta['core_complete'][:,fi]);ri=np.flatnonzero(ref&meta['core_complete'][:,fi]);reference=(C[ri[0],fi],meta['core_denominators'][ri[0],fi]) if len(ri) else None
  summaries.append(run_case(family+'_core',C[use,fi],meta['core_denominators'][use,fi],use,codes,meta,reps=20,reference=reference))
  den=meta['valid_starts'][:,:,fi].sum(axis=1);use=np.flatnonzero(~ref&(den>0));summaries.append(run_case(family+'_all_available',np.asarray(B[:,fi::2][:,use]).T,den[use],use,codes,meta,reps=5))
  for ci,ch in enumerate(meta['chromosomes']):
   den=meta['valid_starts'][:,ci,fi];use=np.flatnonzero(~ref&(den>0))
   if len(use)<100:continue
   summaries.append(run_case(family+'_'+str(ch),CC[use,ci,fi,:],den[use],use,codes[pool[fi]],meta))
 # Diagnose concordance with missingness-sensitive and chromosome-specific analyses.
 comparisons=[];associations=[]
 for fi,family in enumerate(['HSat2','HSat3']):
  with np.load(OUT/'clustering'/(family+'_core')/'model.npz') as z:primary={int(h):int(g) for h,g in zip(z['hap_indices'],z['labels'])};pcs=z['scores'];hi=z['hap_indices']
  for summary in summaries:
   if not summary['case'].startswith(family+'_') or summary['case']==family+'_core':continue
   if not (OUT/'clustering'/summary['case']/'model.npz').exists():comparisons.append([family,summary['case'],0,np.nan]);continue
   with np.load(OUT/'clustering'/summary['case']/'model.npz') as z:other={int(h):int(g) for h,g in zip(z['hap_indices'],z['labels'])}
   common=sorted(set(primary)&set(other));score=ari(np.array([primary[i] for i in common]),np.array([other[i] for i in common])) if len(set(primary.values()))>1 and len(set(other.values()))>1 else np.nan
   comparisons.append([family,summary['case'],len(common),score])
  for variable,values in [('available_chromosomes',meta['available'][hi].sum(axis=1)),('log_core_valid_starts',np.log1p(meta['core_denominators'][hi,fi]))]:
   # Pearson correlation is descriptive, not a significance test of independent subjects.
   for pc in range(min(5,pcs.shape[1])):associations.append([family,pc+1,variable,float(np.corrcoef(pcs[:,pc],values)[0,1]) if np.std(values)>0 else np.nan])
 tsv(OUT/'cluster_concordance.tsv',['family','comparison','overlap_haps','adjusted_rand_index'],comparisons);tsv(OUT/'coverage_associations.tsv',['family','PC','variable','Pearson_r'],associations)
 assignments={}
 for family in ['HSat2','HSat3']:
  with (OUT/'clustering'/(family+'_core')/'assignments.tsv').open() as f:assignments[family]={r['hap']:r for r in csv.DictReader(f,delimiter='\t')}
 tsv(OUT/'hap_groups.tsv',['hap','is_reference','HSat2_group','HSat2_status','HSat3_group','HSat3_status'],([str(h),int(ref[i])]+[value for family in ['HSat2','HSat3'] for value in ([assignments[family][str(h)]['group'],assignments[family][str(h)]['status']] if str(h) in assignments[family] else ['NA','reference_projection_only' if ref[i] else 'incomplete_core_regions'])] for i,h in enumerate(meta['haplotypes'])))
 (OUT/'clustering_summary.json').write_text(json.dumps(summaries,indent=2)+'\n');print('CLUSTERING_COMPLETE',flush=True)
if __name__=='__main__':cluster()
