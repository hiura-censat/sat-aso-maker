#!/usr/bin/env python3
"""Finish scoring, plots and validation when all per-hap counts are complete."""
import json,subprocess,sys,time,hashlib
from pathlib import Path
from step1_pipeline import ROOT,OUT,MANIFEST
while True:
 completed=0
 for row in MANIFEST:
  path=OUT/'counts'/(row['hap']+'.json')
  try:
   meta=json.loads(path.read_text())
   if meta.get('canonical_orientation_agreement') and (OUT/'counts'/(row['hap']+'.npz')).exists():completed+=1
  except (FileNotFoundError,json.JSONDecodeError):pass
 if completed==len(MANIFEST):break
 print(f'Waiting for counts: {completed}/{len(MANIFEST)}',flush=True)
 time.sleep(30)
commands=[([sys.executable,str(ROOT/'scripts/step1_pipeline.py'),'--phase','score'],'scoring.log'),([sys.executable,str(ROOT/'scripts/plot_step1_data.py')],'plot_data.log'),(['Rscript',str(ROOT/'scripts/plot_step1.R')],'plot.log'),([sys.executable,str(ROOT/'scripts/verify_step1.py')],'validation.log')]
for cmd,name in commands:
 print('Running',name,flush=True)
 with (OUT/name).open('w') as log:subprocess.run(cmd,cwd=ROOT,stdout=log,stderr=subprocess.STDOUT,check=True)
validation=json.loads((OUT/'validation_summary.json').read_text());assert validation['status']=='PASS'
validation['finished_at']=time.strftime('%Y-%m-%dT%H:%M:%S%z')
validation['script_sha256']={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in (ROOT/'scripts').glob('*step1*') if p.is_file()}
(OUT/'COMPLETE.json').write_text(json.dumps(validation,indent=2)+'\n')
print('STEP1_COMPLETE',flush=True)
