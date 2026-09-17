#!/usr/bin/env python3
"""Run STEP4 aggregation, plotting, independent audits, and completion marker."""
import csv, gzip, hashlib, json, os, subprocess, sys, time
from pathlib import Path
import numpy as np
from step3_prepare import ROOT
from step2_pipeline import MANIFEST
from step3_neighbors import neighbors, seq, revcode

S3=ROOT/'step3'; S4=ROOT/'step4'; M=S4/'mismatch'; PY=os.environ.get('SAT_ASO_PYTHON',sys.executable)
def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(path):
    with (gzip.open(path,'rt') if str(path).endswith('.gz') else open(path)) as f:return list(csv.DictReader(f,delimiter='\t'))

def main():
    source_sha=sha(S3/'COMPLETE.json')
    subprocess.run([PY,str(ROOT/'scripts/step4_summarize.py')],cwd=ROOT,check=True)
    queries=read(M/'queries.tsv'); new=read(M/'new_queries.tsv'); assert len(queries)>0 and len(new)<=len(queries)
    assert len({r['kmer_id'] for r in queries})==len(queries)
    ranked=read(S4/'ranked/mismatch_evaluated_ranked.tsv'); assert len(ranked)==len(queries) and {r['kmer_id'] for r in ranked}=={r['kmer_id'] for r in queries}
    assert all(r['canonical_target_5to3']==r['aso_5to3_if_reverse_complement_target_is_transcribed'] for r in ranked)
    assert all(r['reverse_complement_target_5to3']==r['aso_5to3_if_canonical_target_is_transcribed'] for r in ranked)
    exact=read(S4/'ranked/all_exact_ranked.tsv.gz'); assert len(exact)>0 and len({r['kmer_id'] for r in exact})==len(exact)
    # Independently reconstruct all radius neighborhoods and 12 evenly spaced saved new-query spectra.
    with np.load(M/'variant_mapping.npz') as z: records=z['records']; vc=z['variant_codes']; qc=z['query_codes']
    assert len(qc)==len(new) and len(records)==len({(int(r[0]),int(r[1])) for r in records})
    assert len(queries)==len(new)+sum(r['query_origin']=='STEP3' for r in queries)
    for qi,q in enumerate(qc):
        rr=records[records[:,1]==qi]; assert set(map(int,vc[rr[:,0]]))==neighbors(int(q))
        for vi,_,d,m1,m2 in rr:
            a=seq(q); b=seq(vc[vi]); rb=seq(revcode(int(vc[vi])))
            assert min(sum(x!=y for x,y in zip(a,b)),sum(x!=y for x,y in zip(a,rb)))==d
            assert int(m1).bit_count()==int(m2).bit_count()==d
    audit_haps=min(12,len(MANIFEST))
    for hi in np.linspace(0,len(MANIFEST)-1,audit_haps,dtype=int):
        h=MANIFEST[hi]['hap']
        with np.load(M/'variants'/(h+'.npz')) as z: raw=z['counts']
        ag=np.zeros((len(qc)*3,raw.shape[1],4),dtype=np.uint64);np.add.at(ag,records[:,1]*3+records[:,2],raw[records[:,0]])
        with np.load(M/'counts'/(h+'.npz')) as z: np.testing.assert_array_equal(ag.reshape(len(qc),3,raw.shape[1],4).transpose(0,2,3,1),z['counts'])
    for r in MANIFEST:
        qcj=json.loads((M/'counts'/(r['hap']+'.json')).read_text())
        assert qcj['zero_mismatch_matches_step1'] and qcj['all_valid_windows_match_step0'] and qcj['mismatch_position_accounting']
    with (S4/'logs/plots.log').open('w') as log: subprocess.run(['Rscript',str(ROOT/'scripts/plot_step4.R')],cwd=ROOT,stdout=log,stderr=subprocess.STDOUT,check=True)
    plots=sorted(S4.glob('plots/*.png')); assert len(plots)==3 and all(p.stat().st_size>10000 for p in plots)
    mm=json.loads((M/'summary.json').read_text()); assert mm['status']=='PASS' and mm['sequence_sets']==len(MANIFEST)
    validation={'status':'PASS','source_step3_complete_sha256':source_sha,'source_step3_unchanged':sha(S3/'COMPLETE.json')==source_sha,'integrated_exact_candidates':len(exact),'mismatch_queries':len(queries),'new_mismatch_queries':len(new),'new_query_radius2_neighborhoods_exhaustive_and_deduplicated':True,'variant_spectrum_serialization_audit_haps':audit_haps,'all_new_query_counts_match_step1_at_0mm':True,'all_valid_windows_match_step0':True,'mismatch_position_accounting_verified':True,'orientation_fields_verified':True,'plots_verified':len(plots)}
    (S4/'validation_summary.json').write_text(json.dumps(validation,indent=2)+'\n')
    nonreference_haps=sum(r['sample']!='CHM13' for r in MANIFEST)
    summary={'status':'complete','integrated_exact_candidates':len(exact),'mismatch_evaluated_candidates':len(queries),'sequence_sets':len(MANIFEST),'nonreference_haps':nonreference_haps,'strict_pan_pass_counts_by_radius':mm['strict_pan_pass_counts_by_radius'],'category_pass_counts_by_radius':mm['category_pass_counts_by_radius'],'robustness_tiers':mm['robustness_tiers'],'scope':'HSat2/HSat3 target discovery using existing annotations on eligible T2T chromosomes','interpretation':'Ranked targets for experimental ASO design; sequence-only quality proxies and genomic occurrence are not knockdown efficacy.'}
    (S4/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    validation['finished_at']=time.strftime('%Y-%m-%dT%H:%M:%S%z'); validation['script_sha256']={p.name:sha(p) for p in (ROOT/'scripts').glob('*step4*') if p.is_file()}
    (S4/'COMPLETE.json').write_text(json.dumps(validation,indent=2)+'\n'); print(json.dumps(summary,indent=2))
if __name__=='__main__':main()
