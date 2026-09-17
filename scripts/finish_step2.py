#!/usr/bin/env python3
import json,subprocess,sys,time
from step2_pipeline import OUT,ROOT,MANIFEST,sha
while True:
 done=sum((OUT/'counts'/(r['hap']+'.json')).exists() and (OUT/'counts'/(r['hap']+'.npz')).exists() for r in MANIFEST)
 if done==len(MANIFEST):break
 print('Waiting for counts',done,len(MANIFEST),flush=True);time.sleep(30)
for command,log in [([sys.executable,'scripts/summarize_step2.py'],'summary.log'),(['Rscript','scripts/plot_step2.R'],'plots.log'),([sys.executable,'scripts/verify_step2.py'],'validation.log')]:
 print('Running',log,flush=True)
 with (OUT/log).open('w') as f:subprocess.run(command,cwd=ROOT,stdout=f,stderr=subprocess.STDOUT,check=True)
result=json.loads((OUT/'validation_summary.json').read_text());assert result['status']=='PASS'
result['finished_at']=time.strftime('%Y-%m-%dT%H:%M:%S%z');result['script_sha256']={p.name:sha(p) for p in (ROOT/'scripts').glob('*step2*') if p.is_file()}
(OUT/'COMPLETE.json').write_text(json.dumps(result,indent=2)+'\n');print('STEP2_COMPLETE',flush=True)
