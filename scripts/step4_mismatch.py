#!/usr/bin/env python3
import argparse
from concurrent.futures import ThreadPoolExecutor,as_completed
import step3_mismatch as engine
from step3_prepare import ROOT
from step2_pipeline import MANIFEST
engine.M=ROOT/'step4/mismatch'
def main(workers,pilot):
 shared=engine.setup();jobs=MANIFEST[:pilot] if pilot else MANIFEST
 with ThreadPoolExecutor(max_workers=workers) as pool:
  for i,f in enumerate(as_completed([pool.submit(engine.count_hap,r,shared) for r in jobs]),1):print('mismatch',i,len(jobs),*f.result(),flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--workers',type=int,default=8);p.add_argument('--pilot',type=int,default=0);a=p.parse_args();main(a.workers,a.pilot)
