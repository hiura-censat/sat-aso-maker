import json,sys
import numpy as np
import step1_pipeline as p
original=p.OUT;work=original/'prefix_validation';work.mkdir(exist_ok=True);(work/'logs').mkdir(exist_ok=True)
row=next(r for r in p.MANIFEST if r['hap']=='HG00096_hap1');codes=np.fromfile(original/'preliminary_candidates.u32',dtype='<u4')
with np.load(original/'counts/HG00096_hap1.npz') as z:expected=z['counts'][:,3,:]
p.OUT=work;rawfile=work/'background.bin.gz'
qc=p.count_region(row,'background',rawfile,canonical=False,seed=original/'preliminary_candidates.both_strands.fa')
raw=p.readbin(rawfile);rev=p.revcodes(raw['code']);can=np.minimum(raw['code'],rev);indices=np.searchsorted(codes,can)
assert np.all(indices<len(codes)) and np.array_equal(codes[indices],can)
actual=np.zeros_like(expected);actual[indices,(raw['code']>rev).astype(int)]=raw['count']
assert np.array_equal(actual,expected)
qc.update(hap=row['hap'],all_background_candidate_orientation_counts_identical_to_stock=True)
(original/'prefix_validation.json').write_text(json.dumps(qc,indent=2)+'\n');rawfile.unlink();print(json.dumps(qc),flush=True)
