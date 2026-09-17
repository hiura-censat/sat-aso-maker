import collections,random,subprocess,tempfile,unittest
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
def rc(s):return s.translate(str.maketrans('ACGT','TGCA'))[::-1]
def encode(s):
 n=0
 for b in s:n=(n<<2)|'ACGT'.index(b)
 return n
class ExactCounter(unittest.TestCase):
 def test_boundaries_strands_ambiguous_lowercase_chromosomes(self):
  rng=random.Random(81)
  records=[('chr1','a'*70000),('chr2','T'*70001),('chr1','ACGT'*10),('chr2','AAA'),('chr2','A'*13),('chr1',''.join(rng.choices('ACGTNacgt',k=7000))),('chr2','GATTACA'*40)]
  full=collections.Counter();valid=collections.Counter()
  for chrom,s in records:
   s=s.upper()
   for i in range(len(s)-15):
    k=s[i:i+16]
    if set(k)<=set('ACGT'):full[(min(k,rc(k)),chrom,int(k>rc(k)))]+=1;valid[chrom]+=1
  keys=sorted({x[0] for x in full});keys=sorted(set(keys[::2])|{'A'*16,'C'*16,'ACGT'*4})
  with tempfile.TemporaryDirectory() as tmp:
   p=Path(tmp);np.array([encode(k) for k in keys],dtype='<u4').tofile(p/'codes');(p/'chroms').write_text('chr1\nchr2\nchrX\n')
   fa=''.join(f'>{c}:{i+1}-{i+len(s)}\n'+'\n'.join(s[j:j+61] for j in range(0,len(s),61))+'\n' for i,(c,s) in enumerate(records))
   subprocess.run([str(ROOT/'scripts/step2_count'),str(p/'codes'),str(p/'chroms'),str(p/'out'),str(p/'stats')],input=fa,text=True,check=True)
   actual=np.fromfile(p/'out',dtype='<u8').reshape(len(keys),3,2)
   for i,k in enumerate(keys):
    for j,c in enumerate(['chr1','chr2','chrX']):
     for o in range(2):self.assertEqual(int(actual[i,j,o]),full[(k,c,o)])
   import csv
   with (p/'stats').open() as f:
    for row in csv.DictReader(f,delimiter='\t'):self.assertEqual(int(row['valid_16mer_starts']),valid[row['chromosome']])
if __name__=='__main__':unittest.main()
