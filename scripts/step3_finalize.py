#!/usr/bin/env python3
import csv,hashlib,json,subprocess,time
from pathlib import Path
import numpy as np
from step3_prepare import ROOT,OUT,S2
from step3_neighbors import neighbors,seq,revcode
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def finalize():
 config=json.loads((OUT/'config.json').read_text());assert sha(S2/'COMPLETE.json')==config['source_complete_sha256']
 for path,old in json.loads((ROOT/'step0/input_metadata.json').read_text()).items():
  p=Path(path)
  if old is None:assert not p.exists()
  else:assert {'size':p.stat().st_size,'mtime_ns':p.stat().st_mtime_ns}==old,path
 with np.load(OUT/'data/metadata.npz') as z:meta={k:z[k] for k in z.files}
 summaries=json.loads((OUT/'clustering_summary.json').read_text())
 chromosome_counts=np.load(OUT/'data/chromosome_feature_counts.npy',mmap_mode='r')
 feature_pool=np.load(OUT/'data/chromosome_feature_indices.npy')
 for case in summaries:
  if case['haps']<100 or case.get('status')=='insufficient_variable_features_in_analyzed_pool':continue
  fi=['HSat2','HSat3'].index(case['family']);ci=list(meta['chromosomes']).index(case['chromosome'])
  folder=OUT/'clustering'/case['case']
  with np.load(folder/'model.npz') as z:
   hi=z['hap_indices'];codes=z['feature_codes'];labels=z['labels'];features=z['features']
   pool_codes=meta['codes'][feature_pool[ci,fi]]
   ix=np.searchsorted(np.sort(pool_codes),codes)
   assert np.array_equal(np.sort(pool_codes)[ix],codes)
   source_index=np.array([np.flatnonzero(pool_codes==code)[0] for code in codes])
   assert not meta['is_reference'][hi].any() and np.all(meta['valid_starts'][hi,ci,fi]>0)
   expected=np.log1p(np.asarray(chromosome_counts[hi,ci,fi,:])[:,source_index]/meta['valid_starts'][hi,ci,fi,None]*1e6)
   assert np.allclose(features,expected,rtol=1e-12,atol=1e-12)
   if features.shape[1]>1:
    corr=np.corrcoef(features.T);np.fill_diagonal(corr,0);assert np.nanmax(np.abs(corr))<.950001
   assert len(np.unique(labels))==case['selected_k'] and len(hi)==case['haps']
   if case['supported_partition']:
    assert min(np.bincount(labels)[1:])>=10
    assert float(case['selection']['silhouette'])>=.25 and float(case['selection']['median_subsample_ARI'])>=.75 and float(case['selection']['median_worst_cluster_Jaccard'])>=.75
 groups=list(csv.DictReader((OUT/'hap_chromosome_groups.tsv').open(),delimiter='\t'))
 assert len(groups)==574*2*len(meta['chromosomes'])
 lookup={(r['hap'],r['family'],r['chromosome']):r for r in groups}
 assert len(lookup)==len(groups)
 for hi,hap in enumerate(meta['haplotypes']):
  for fi,fam in enumerate(['HSat2','HSat3']):
   for ci,ch in enumerate(meta['chromosomes']):
    r=lookup[(str(hap),fam,str(ch))]
    if meta['is_reference'][hi] or meta['valid_starts'][hi,ci,fi]==0:assert r['group']=='NA'
    if r['group']!='NA':assert r['measurable']=='1' and r['status']=='stable_partition'
 with np.load(OUT/'mismatch/variant_mapping.npz') as z:records=z['records'];vc=z['variant_codes'];qc=z['query_codes']
 assert len(records)==len({(int(r[0]),int(r[1])) for r in records})
 for qi,q in enumerate(qc):
  r=records[records[:,1]==qi];assert set(map(int,vc[r[:,0]]))==neighbors(int(q))
  for vi,_,d,m1,m2 in r:
   qseq=seq(q);v=seq(vc[vi]);rv=seq(revcode(int(vc[vi])));a=sum(x!=y for x,y in zip(qseq,v));b=sum(x!=y for x,y in zip(qseq,rv));assert min(a,b)==d;assert int(m1).bit_count()==int(m2).bit_count()==d
 # Audit saved variants against spectra for 12 evenly spaced haps, independently of serialization.
 for hi in np.linspace(0,573,12,dtype=int):
  h=str(meta['haplotypes'][hi]);path=OUT/'mismatch'
  with np.load(path/'variants'/(h+'.npz')) as z:c=z['counts']
  aggregate=np.zeros((len(qc)*3,c.shape[1],4),dtype=np.uint64);np.add.at(aggregate,records[:,1]*3+records[:,2],c[records[:,0]])
  with np.load(path/'counts'/(h+'.npz')) as z:assert np.array_equal(aggregate.reshape(len(qc),3,c.shape[1],4).transpose(0,2,3,1),z['counts'])
 for h in meta['haplotypes']:
  m=json.loads((OUT/'mismatch/counts'/(str(h)+'.json')).read_text());assert m['zero_mismatch_matches_step1'] and m['all_valid_windows_match_step0'] and m['mismatch_position_accounting']
 mm=json.loads((OUT/'mismatch/summary.json').read_text());assert mm['status']=='PASS' and mm['sequence_sets']==574
 with (OUT/'logs/plots.log').open('w') as log:subprocess.run(['Rscript',str(ROOT/'scripts/plot_step3.R')],cwd=ROOT,stdout=log,stderr=subprocess.STDOUT,check=True)
 assert len(list((OUT/'plots').glob('*.png')))>=2
 validation={'status':'PASS','source_metadata_unchanged':True,'chromosome_clusters_reconstructed_from_measurable_haps':True,'reference_excluded_from_fitting':True,'missing_chromosome_haps_unassigned':True,'representative_correlations_below_0.95':True,'cluster_acceptance_criteria_verified':True,'all_radius_2_neighborhoods_exhaustive_and_deduplicated':True,'variant_spectrum_serialization_audit_haps':12,'all_574_hap_exact_counts_match_step1':True,'mismatch_position_accounting_verified':True}
 (OUT/'validation_summary.json').write_text(json.dumps(validation,indent=2)+'\n')
 summary={'status':'complete','chromosome_family_cases':len(summaries),'supported_chromosome_partitions':sum(x['supported_partition'] for x in summaries),'mismatch_queries':mm['queries'],'mismatch_sequence_sets':574,'strict_pan_pass_counts_by_radius':mm['strict_pan_pass_counts_by_radius'],'interpretation':'Independent family-by-chromosome groups; no whole-genome hap group is defined.'}
 (OUT/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');validation['finished_at']=time.strftime('%Y-%m-%dT%H:%M:%S%z');validation['script_sha256']={p.name:sha(p) for p in (ROOT/'scripts').glob('*step3*') if p.is_file()};(OUT/'COMPLETE.json').write_text(json.dumps(validation,indent=2)+'\n');print(json.dumps(summary,indent=2),flush=True)
if __name__=='__main__':finalize()
