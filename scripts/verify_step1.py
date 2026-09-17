#!/usr/bin/env python3
import csv,hashlib,json,sys
from pathlib import Path
import numpy as np
from step1_pipeline import ROOT,OUT,S0,CONFIG,MANIFEST,REGIONS,STATS,fingerprint,revcodes

def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
for path,expected in json.loads((OUT/'input_checksums.json').read_text()).items():assert sha(ROOT/path)==expected,path
for path,old in json.loads((S0/'input_metadata.json').read_text()).items():
 p=Path(path)
 if old is None:assert not p.exists(),path
 else:assert {'size':p.stat().st_size,'mtime_ns':p.stat().st_mtime_ns}==old,path
codes=np.fromfile(OUT/'preliminary_candidates.u32',dtype='<u4');assert len(codes)>0
assert np.all(np.diff(codes.astype(np.uint64))>0)
assert np.all(codes<=revcodes(codes))
seedhash=sha(OUT/'preliminary_candidates.u32')
for r in MANIFEST:
 h=r['hap'];d=json.loads((OUT/'discovery'/(h+'.bin.json')).read_text())
 assert d['fingerprint']==fingerprint(r,'HSat23')
 assert d['count_sum']==STATS[(h,'HSat23')]
 qc=json.loads((OUT/'counts'/(h+'.json')).read_text())
 expected=hashlib.sha256(''.join(fingerprint(r,f,seedhash) for f in REGIONS).encode()).hexdigest()
 assert qc['fingerprint']==expected,h
 assert qc['canonical_orientation_agreement']
 for f in REGIONS:
  row=qc['regions'][f];assert row['expected_valid_starts']==STATS[(h,f)]<2**32
  assert row['candidate_count_sum']<=row['count_sum']<=row['expected_valid_starts']
  if f!='background':assert row['count_sum']==row['expected_valid_starts']
 assert (OUT/'enrichment_by_hap'/(h+'.npz')).exists()
with np.load(OUT/'pooled_counts_and_scores.npz') as z:
 assert np.array_equal(z['codes'],codes)
 counts=z['counts'];t=counts[:,2,:].sum(axis=1);b=counts[:,3,:].sum(axis=1);den=z['valid_starts']
 expected_den=np.array([sum(STATS[(r['hap'],f)] for r in MANIFEST if r['sample']!='CHM13') for f in REGIONS],dtype=np.uint64)
 assert np.array_equal(den,expected_den)
 expected=np.log2((t/den[2]+1e-9)/(b/den[3]+1e-9));assert np.allclose(z['log2_enrichment'],expected,rtol=0,atol=1e-12)
 pooled=(t>=100)&(expected>=2);local=z['qualifying_haps']>0;selected=pooled|local
 assert np.array_equal(z['selected'],selected)
 assert np.array_equal(np.fromfile(OUT/'candidate_codes.u32',dtype='<u4'),codes[selected])
 assert np.all(z['prevalence_haps']<=573)
 assert np.all(z['qualifying_haps']<=z['prevalence_haps'])
 expected_count=int(selected.sum())
with (OUT/'candidate_kmers.tsv').open() as f:
 seen=set()
 for row in csv.DictReader(f,delimiter='\t'):
  k=row['canonical_kmer'];assert len(k)==16 and set(k)<=set('ACGT')
  reverse=k.translate(str.maketrans('ACGT','TGCA'))[::-1];assert reverse==row['reverse_complement'] and k<=reverse
  code=0
  for ch in k:code=(code<<2)|'ACGT'.index(ch)
  assert row['kmer_id']==f'k16_{code:08x}' and code not in seen;seen.add(code)
 assert len(seen)==expected_count
 assert seen==set(map(int,codes[selected]))
summary={'status':'PASS','sequence_sets':len(MANIFEST),'nonreference_haps':573,'preliminary_candidates':len(codes),'selected_candidates':expected_count,'all_target_counts_match_step0':True,'canonical_orientation_agreement':True,'source_metadata_unchanged':True,'selection_and_scores_verified':True,'background_counts_are_candidate_restricted':True,'python':sys.version,'jellyfish':'2.3.1'}
(OUT/'validation_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps(summary,indent=2))
