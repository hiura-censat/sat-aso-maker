#!/usr/bin/env python3
import os,subprocess,sys
from step3_prepare import ROOT,OUT
env=os.environ.copy();env['OPENBLAS_NUM_THREADS']='2';env['OMP_NUM_THREADS']='2'
stages=[('step3_cluster.py',[],'clustering.log'),('step3_candidates.py',[],'candidates.log'),('step3_neighbors.py',[],'neighbors.log'),('step3_mismatch.py',['--workers','2','--pilot','2'],'mismatch_pilot.log'),('step3_mismatch.py',['--workers',os.environ.get('SAT_ASO_STEP3_WORKERS','8')],'mismatch.log'),('step3_summarize_mismatch.py',[],'mismatch_summary.log'),('step3_finalize.py',[],'finalize.log')]
for script,args,log in stages:
 print('Running',script,flush=True)
 with (OUT/'logs'/log).open('w') as f:subprocess.run([sys.executable,str(ROOT/'scripts'/script)]+args,cwd=ROOT,env=env,stdout=f,stderr=subprocess.STDOUT,check=True)
print('STEP3_COMPLETE',flush=True)
