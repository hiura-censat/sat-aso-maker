#!/usr/bin/env python3
import csv,json
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'step1';P=OUT/'plots';P.mkdir(exist_ok=True)
def write(name,header,rows):
 with (P/name).open('w') as f:
  w=csv.writer(f,delimiter='\t');w.writerow(header);w.writerows(rows)
summary=json.loads((OUT/'summary.json').read_text());discovery=json.loads((OUT/'discovery_summary.json').read_text())
write('summary.tsv',['metric','value'],list(summary.items())+list(discovery.items()))
with np.load(OUT/'pooled_counts_and_scores.npz') as z:
 e=z['log2_enrichment'];selected=z['selected'];codes=z['codes'];target=z['target_density'];background=z['background_density'];prev=z['prevalence_haps'];counts=z['counts']
 hist,edges=np.histogram(e,bins=np.linspace(np.floor(e.min()),np.ceil(e.max())+1e-6,101))
 write('enrichment_histogram.tsv',['midpoint','count'],zip((edges[:-1]+edges[1:])/2,hist))
 rng=np.random.default_rng(42);indices=rng.choice(len(codes),min(50000,len(codes)),replace=False)
 write('density_scatter.tsv',['log10_background_plus_epsilon','log10_target_plus_epsilon','selected'],zip(np.log10(background[indices]+1e-9),np.log10(target[indices]+1e-9),selected[indices].astype(int)))
 x=codes[selected].copy();gc=np.zeros(len(x),dtype=int)
 for _ in range(16):b=x&3;gc+=(b==1)|(b==2);x>>=2
 write('gc_histogram.tsv',['GC_percent','count'],zip(np.arange(17)*100/16,np.bincount(gc,minlength=17)))
 h,edges=np.histogram(prev[selected]/summary['nonreference_haps'],bins=np.linspace(0,1,21));write('prevalence_histogram.tsv',['midpoint','count'],zip((edges[:-1]+edges[1:])/2,h))
 forward=counts[selected,2,0].astype(float);total=counts[selected,2,:].sum(axis=1);frac=forward/total
 h,edges=np.histogram(frac,bins=np.linspace(0,1,21));write('orientation_histogram.tsv',['midpoint','count'],zip((edges[:-1]+edges[1:])/2,h))
