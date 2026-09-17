import collections,subprocess,tempfile,unittest
from pathlib import Path
import numpy as np
from step3_neighbors import neighbors,encode,revcode,seq,distance_masks
from step3_cluster import ari
ROOT=Path(__file__).resolve().parents[1]
class Step3Tests(unittest.TestCase):
 def test_adjusted_rand(self):
  self.assertAlmostEqual(ari(np.array([1,1,2,2]),np.array([9,9,3,3])),1)
  self.assertAlmostEqual(ari(np.array([1,1,2,2]),np.array([1,2,1,2])),-.5)
  self.assertEqual(ari(np.ones(4),np.ones(4)),1)
 def test_neighbors_counter_positions(self):
  queries=['ACGTTGCAATGCCGTA','ACGT'*4,'A'*16];codes=[min(encode(s),revcode(encode(s))) for s in queries];queries=[seq(c) for c in codes]
  variants=np.array(sorted(set.union(*(neighbors(c) for c in codes))),dtype='<u4')
  self.assertEqual(len(neighbors(codes[0])),1129)
  records=['a'*70000,'T'*60,'ACGTTGCAATGCCGTA','ACGTTGCAATGCCGTT','ACGTTGCAATGCCGTC','ACGT'*20,'NNNNACGTACGTACGTACGTNNNN','ACGTTGCA','ATGCCGTA']
  want=np.zeros((len(codes),3),dtype=np.uint64)
  for s in records:
   s=s.upper()
   for i in range(len(s)-15):
    k=s[i:i+16]
    if set(k)-set('ACGT'):continue
    kr=seq(revcode(encode(k)))
    for qi,q in enumerate(queries):
     d=min(sum(a!=b for a,b in zip(q,k)),sum(a!=b for a,b in zip(q,kr)))
     if d<=2:want[qi,d]+=1
  with tempfile.TemporaryDirectory() as tmp:
   p=Path(tmp);variants.tofile(p/'codes');(p/'chr').write_text('chr1\n');fa=''.join(f'>chr1:{i+1}-{i+len(s)}\n{s}\n' for i,s in enumerate(records))
   subprocess.run([str(ROOT/'scripts/step3_count'),str(p/'codes'),str(p/'chr'),str(p/'counts'),str(p/'stats')],input=fa,text=True,check=True)
   counted=np.fromfile(p/'counts',dtype='<u8').reshape(len(variants),1,2).sum(axis=(1,2));got=np.zeros_like(want)
   for qi,q in enumerate(codes):
    for v in neighbors(q):
     d,m1,m2=distance_masks(q,v);assert m1.bit_count()==m2.bit_count()==d;got[qi,d]+=counted[np.searchsorted(variants,v)]
   np.testing.assert_array_equal(got,want)
if __name__=='__main__':unittest.main()
