import collections,subprocess,tempfile,unittest
from pathlib import Path
import numpy as np
from step3_neighbors import neighbors,encode,revcode,seq,distance_masks
from step3_cluster import ari,measurable_haps
import step3_cluster as cluster_module
from step3_prepare import choose_feature_pool
ROOT=Path(__file__).resolve().parents[1]
class Step3Tests(unittest.TestCase):
 def test_feature_pool_is_chromosome_independent(self):
  counts=np.zeros((10,2,2),dtype=np.uint64)
  counts[:,0,0]=np.arange(10)
  original=choose_feature_pool(counts,width=6)[0,0].copy()
  counts[:,1,0]=np.arange(10)[::-1]*100
  np.testing.assert_array_equal(choose_feature_pool(counts,width=6)[0,0],original)
 def test_chromosome_eligibility_ignores_other_chromosomes(self):
  den=np.array([[[100,0],[0,100]],[[100,0],[100,100]],[[0,0],[100,100]],[[100,0],[100,100]]])
  ref=np.array([False,False,False,True])
  np.testing.assert_array_equal(measurable_haps(den,ref,0,0),[0,1])
  np.testing.assert_array_equal(measurable_haps(den,ref,1,0),[1,2])
  den[0,1,0]=1000
  np.testing.assert_array_equal(measurable_haps(den,ref,0,0),[0,1])
 def test_cluster_assigns_only_measurable_chromosome_haps(self):
  with tempfile.TemporaryDirectory() as td:
   old=cluster_module.OUT;cluster_module.OUT=Path(td)
   try:
    (Path(td)/'data').mkdir();(Path(td)/'clustering').mkdir()
    nh=102;den=np.zeros((nh,2,2),dtype=np.uint32);den[:100,0,0]=1000;den[100,1,0]=1000;den[101,0,0]=1000
    ref=np.zeros(nh,dtype=bool);ref[101]=True
    np.savez_compressed(Path(td)/'data/metadata.npz',codes=np.arange(5000,dtype=np.uint32),haplotypes=np.array([f'h{i}' for i in range(nh)]),is_reference=ref,chromosomes=np.array(['chr1','chr2']),valid_starts=den,available=den[:,:,0]>0)
    np.save(Path(td)/'data/chromosome_feature_indices.npy',np.tile(np.arange(5000,dtype=np.int32),(2,2,1)))
    counts=np.lib.format.open_memmap(Path(td)/'data/chromosome_feature_counts.npy',mode='w+',dtype=np.uint32,shape=(nh,2,2,5000));counts[:]=0
    rng=np.random.default_rng(5);counts[:100,0,0,:20]=rng.poisson(15+np.repeat([0,1],50)[:,None]*np.arange(1,21)[None,:]*3);counts.flush();del counts
    cluster_module.cluster()
    with (Path(td)/'hap_chromosome_groups.tsv').open() as f:rows=list(__import__('csv').DictReader(f,delimiter='\t'))
    lookup={(r['hap'],r['family'],r['chromosome']):r for r in rows}
    self.assertEqual(len(rows),nh*2*2)
    self.assertTrue(lookup['h0','HSat2','chr1']['group'].startswith('G'))
    self.assertEqual(lookup['h100','HSat2','chr1']['status'],'missing_target_region')
    self.assertEqual(lookup['h100','HSat2','chr2']['status'],'insufficient_measurable_haps')
    self.assertEqual(lookup['h101','HSat2','chr1']['group'],'NA')
   finally:cluster_module.OUT=old
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
