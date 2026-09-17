"""Complete remaining sequence scans with the compiled scanner, in reverse order."""
import csv, subprocess
from concurrent.futures import ThreadPoolExecutor,as_completed
from step0_preprocess import ROOT, OUT
with (OUT/'manifest.tsv').open() as f:
    jobs=[r for r in csv.DictReader(f,delimiter='\t') if r['status']=='PASS']
jobs=jobs[len(jobs)//2:]+jobs[:len(jobs)//2]
def run(r):
    hap=r['hap']; dest=OUT/'qc'/f'{hap}.stats.tsv'
    if dest.exists(): return hap,'cached'
    dest.parent.mkdir(exist_ok=True)
    tmp=dest.with_suffix('.extra.tmp')
    subprocess.run([str(ROOT/'scripts/step0_scan'),r['fasta'],str(OUT/'regions'/hap),hap,str(tmp)],check=True)
    tmp.replace(dest)
    return hap,'done'
with ThreadPoolExecutor(max_workers=8) as pool:
    for i,ft in enumerate(as_completed([pool.submit(run,r) for r in jobs]),1):
        print(i,*ft.result(),flush=True)
