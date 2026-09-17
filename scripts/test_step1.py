import collections,subprocess,tempfile,unittest
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];JF=ROOT/'.local/bin/jellyfish'
def rc(s):return s.translate(str.maketrans('ACGT','TGCA'))[::-1]
class JellyfishSemantics(unittest.TestCase):
 def test_exact_canonical_and_selected_strands(self):
  records=['a'*200000,'T'*50,'ACGT'*4,'GATTACA'*10+'N'+'ccgat'*20,'AAAA','TTTTTTTTTTTT']
  expected=collections.Counter()
  for seq in records:
   seq=seq.upper()
   for i in range(len(seq)-15):
    k=seq[i:i+16]
    if set(k)<=set('ACGT'):expected[k]+=1
  canonical=collections.Counter()
  for k,n in expected.items():canonical[min(k,rc(k))]+=n
  with tempfile.TemporaryDirectory() as temp:
   d=Path(temp);fa=d/'input.fa';fa.write_text(''.join(f'>{i}\n{s}\n' for i,s in enumerate(records)))
   seed=['A'*16,'T'*16,'ACGT'*4,'C'*16];sf=d/'seed.fa';sf.write_text(''.join(f'>{i}\n{s}\n' for i,s in enumerate(seed)))
   for flags,want in [([],expected),(['-C'],canonical),(['--if',str(sf)],{k:n for k,n in expected.items() if k in seed})]:
    db=d/'out.jf'
    subprocess.run([str(JF),'count','-m','16','-s','1M','-c','16','--out-counter-len','4','-L','1','-o',str(db)]+flags+[str(fa)],check=True)
    text=subprocess.check_output([str(JF),'dump','-c',str(db)],text=True)
    got={a:int(b) for a,b in (l.split() for l in text.splitlines()) if int(b)>0}
    self.assertEqual(got,dict(want))


class RegionPipeline(unittest.TestCase):
 def test_region_counts_against_bruteforce(self):
  import sys,json
  import numpy as np
  sys.path.insert(0,str(ROOT/'scripts'))
  import step1_pipeline as p
  old_out,old_s0=p.OUT,p.S0
  with tempfile.TemporaryDirectory() as temp:
   d=Path(temp);p.OUT=d/'out';p.S0=d/'s0';hap='fixture'
   try:
    for sub in ['discovery','counts','work','logs']: (p.OUT/sub).mkdir(parents=True)
    beddir=p.S0/'regions'/hap;beddir.mkdir(parents=True)
    seq='A'*70000+'T'*40+'ACGT'*30+'N'+'GATTACA'*25+'ccgat'*20
    fa=d/'input.fa';fa.write_text('>chr1\n'+seq+'\n')
    subprocess.run([str(p.SAM),'faidx',str(fa)],check=True)
    ranges={'HSat2':[(0,1000),(70000,70056)],'HSat3':[(1000,2000),(70056,70160)],'HSat23':[(0,2000),(70000,70160)],'background':[(2000,70000),(70160,len(seq))]}
    wants={}
    for family,ivs in ranges.items():
     (beddir/(family+'.bed')).write_text(''.join(f'chr1\t{s}\t{e}\t{family}\n' for s,e in ivs))
     counter=collections.Counter()
     for s,e in ivs:
      text=seq[s:e].upper()
      for i in range(len(text)-15):
       k=text[i:i+16]
       if set(k)<=set('ACGT'):counter[k]+=1
     wants[family]=counter;p.STATS[(hap,family)]=sum(counter.values())
    row={'hap':hap,'fasta':str(fa),'sample':'fixture'}
    p.discovery(row)
    candidates=sorted({min(k,rc(k)) for k in wants['HSat23']})
    def code(k):
     n=0
     for ch in k:n=(n<<2)|'ACGT'.index(ch)
     return n
    codes=np.array([code(k) for k in candidates],dtype=np.uint32)
    seed=d/'seed.fa';both=sorted(set(candidates+[rc(k) for k in candidates]));seed.write_text(''.join(f'>{i}\n{k}\n' for i,k in enumerate(both)))
    p.candidate_counts(row,codes,seed,'fixture')
    with np.load(p.OUT/'counts'/(hap+'.npz')) as z:actual=z['counts']
    for i,k in enumerate(candidates):
     for ri,family in enumerate(p.REGIONS):
      self.assertEqual(int(actual[i,ri,0]),wants[family][k])
      self.assertEqual(int(actual[i,ri,1]),0 if k==rc(k) else wants[family][rc(k)])
   finally:p.OUT,p.S0=old_out,old_s0
if __name__=='__main__':unittest.main()
