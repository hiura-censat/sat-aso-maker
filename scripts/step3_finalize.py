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
 core=np.load(OUT/'data/core_counts.npy',mmap_mode='r');cluster_summaries=[]
 for fi,family in enumerate(['HSat2','HSat3']):
  folder=OUT/'clustering'/(family+'_core');s=json.loads((folder/'summary.json').read_text());cluster_summaries.append(s)
  with np.load(folder/'model.npz') as z:
   hi=z['hap_indices'];codes=z['feature_codes'];labels=z['labels'];features=z['features'];ix=np.searchsorted(meta['codes'],codes)
   assert np.array_equal(meta['codes'][ix],codes);assert not meta['is_reference'][hi].any();assert meta['core_complete'][hi,fi].all()
   expected=np.log1p(np.asarray(core[hi,fi,:])[:,ix]/meta['core_denominators'][hi,fi,None]*1e6);assert np.allclose(features,expected,rtol=1e-12,atol=1e-12)
   corr=np.corrcoef(features.T);np.fill_diagonal(corr,0);assert np.nanmax(np.abs(corr))<.950001
   assert len(np.unique(labels))==s['selected_k'];assert len(hi)==s['haps']
   if s['supported_partition']:
    assert min(np.bincount(labels)[1:])>=10;assert float(s['selection']['silhouette'])>=.25;assert float(s['selection']['median_subsample_ARI'])>=.75;assert float(s['selection']['median_worst_cluster_Jaccard'])>=.75
 groups=list(csv.DictReader((OUT/'hap_groups.tsv').open(),delimiter='\t'));assert len(groups)==574 and len({r['hap'] for r in groups})==574
 for hi,r in enumerate(groups):
  assert r['hap']==meta['haplotypes'][hi]
  for fi,family in enumerate(['HSat2','HSat3']):
   if not meta['core_complete'][hi,fi] or meta['is_reference'][hi]:assert r[family+'_group']=='NA'
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
 assert len(list((OUT/'plots').glob('*.png')))==5
 validation={'status':'PASS','source_metadata_unchanged':True,'primary_clusters_reconstructed_from_core_counts':True,'reference_excluded_from_fitting':True,'missing_haps_unassigned':True,'representative_correlations_below_0.95':True,'cluster_acceptance_criteria_verified':True,'all_radius_2_neighborhoods_exhaustive_and_deduplicated':True,'variant_spectrum_serialization_audit_haps':12,'all_574_hap_exact_counts_match_step1':True,'mismatch_position_accounting_verified':True}
 (OUT/'validation_summary.json').write_text(json.dumps(validation,indent=2)+'\n')
 summary={'status':'complete','HSat2_core_haps':cluster_summaries[0]['haps'],'HSat2_core_groups':cluster_summaries[0]['selected_k'],'HSat3_core_haps':cluster_summaries[1]['haps'],'HSat3_core_groups':cluster_summaries[1]['selected_k'],'mismatch_queries':mm['queries'],'mismatch_sequence_sets':574,'strict_pan_pass_counts_by_radius':mm['strict_pan_pass_counts_by_radius'],'interpretation':'Regional sequence-pattern groups; whole-genome group concordance is low. Group marker selection is exploratory on the discovery data.'}
 (OUT/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');validation['finished_at']=time.strftime('%Y-%m-%dT%H:%M:%S%z');validation['script_sha256']={p.name:sha(p) for p in (ROOT/'scripts').glob('*step3*') if p.is_file()};(OUT/'COMPLETE.json').write_text(json.dumps(validation,indent=2)+'\n');print(json.dumps(summary,indent=2),flush=True)
if __name__=='__main__':finalize()
