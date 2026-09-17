import sys,time,subprocess
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent))
from step1_pipeline import ROOT,OUT,SAM,JF,MANIFEST
r=next(r for r in MANIFEST if r['hap']=='HG00096_hap1')
fa=OUT/'work/bloom_benchmark.fa'
with fa.open('wb') as f:subprocess.run([str(SAM),'faidx','-n','1000000',r['fasta'],'chr11'],stdout=f,check=True)
def run(use):
 db=OUT/'work'/('benchmark_'+str(use)+'.jf')
 cmd=[str(JF),'count','-m','16','-s','8M','-t','2','-c','16','--out-counter-len','4','-L','1','--if',str(OUT/'preliminary_candidates.both_strands.fa'),'-o',str(db)]
 if use:cmd+=['--bc',str(OUT/'preliminary_candidates.bc')]
 cmd+=[str(fa)];t=time.monotonic();subprocess.run(cmd,check=True);elapsed=time.monotonic()-t
 text=subprocess.check_output([str(JF),'dump','-c',str(db)],text=True);counts={a:int(b) for a,b in (l.split() for l in text.splitlines()) if int(b)>0};db.unlink()
 print('bloom',use,'seconds',elapsed,'nonzero_kmers',len(counts),flush=True)
 return counts
with ThreadPoolExecutor(max_workers=2) as pool:a,b=list(pool.map(run,[False,True]))
assert a==b
print('COUNTS_IDENTICAL',flush=True);fa.unlink()
