#!/usr/bin/env python3
"""Classify all 72 CHM13 exact hits and export hap/group/chromosome count views."""
import csv,gzip,json,bisect,hashlib
from collections import defaultdict
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'step4/all72_landscape';OUT.mkdir(parents=True,exist_ok=True)
FINE=ROOT.parent/'VallePrep_v0.0.0/workdir/work1/bed_all/chm13v2.0.censat.bed';BROAD=ROOT.parent/'VallePrep_v0.0.0/reference/chm13v2.0_censat_v2.1.bed'
CHRS=[f'chr{i}' for i in range(1,23)]+['chrX','chrY'];LABELS=['HSat3','HSat2','alphaSat','ct','outside']
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(path):
 with (gzip.open(path,'rt') if str(path).endswith('.gz') else open(path)) as f:return list(csv.DictReader(f,delimiter='\t'))
def bed(path):
 out=defaultdict(lambda:defaultdict(list))
 for line in path.read_text().splitlines():
  if not line or line.startswith('#'):continue
  x=line.split();out[x[0]][x[3]].append((int(x[1]),int(x[2])))
 return {ch:{label:(tuple(s for s,e in sorted(rs)),tuple(e for s,e in sorted(rs))) for label,rs in by.items()} for ch,by in out.items()}
def match(by,ch,s,e,priority):
 if ch not in by:return 'outside'
 for label in priority:
  if label not in by[ch]:continue
  starts,ends=by[ch][label];j=bisect.bisect_right(starts,s)-1
  while j>=0:
   if e<=ends[j]:return label
   j-=1
 return 'outside'
def write(path,rows,fields):
 with (gzip.open(path,'wt') if str(path).endswith('.gz') else open(path,'w')) as f:
  w=csv.DictWriter(f,fieldnames=fields,delimiter='\t');w.writeheader();w.writerows(rows)
def main():
 q=read(ROOT/'step4/mismatch/queries.tsv');ids=[r['kmer_id'] for r in q];fams={r['kmer_id']:r['target_family'] for r in q};assert ids and len(set(ids))==len(ids)
 fine=bed(FINE);broad=bed(BROAD);assert set(fine)==set(CHRS) and set(broad)==set(CHRS)
 for fam in ['HSat2','HSat3']:
  reference={(ch,s,e) for ch,by in fine.items() for s,e in zip(*by.get(fam,((),())))}
  step0={(x[0],int(x[1]),int(x[2])) for line in (ROOT/'step0/regions/chm13v2.0'/f'{fam}.bed').read_text().splitlines() if (x:=line.split())}
  assert reference==step0,(fam,len(reference),len(step0))
 counts=defaultdict(int);bins=defaultdict(int);raw=OUT/'raw/chm13_all72_exact_hits.tsv.gz';dest=OUT/'chm13_all72_exact_hits_with_censat.tsv.gz';n=0
 fields=['kmer_id','chromosome','start_0based','end_0based','strand_of_canonical_target','target_family','inside_broad_CenSat','fine_censat_family','outside_broad_CenSat','outside_fine_annotation']
 with gzip.open(raw,'rt') as inp,gzip.open(dest,'wt') as out:
  r=csv.DictReader(inp,delimiter='\t');w=csv.DictWriter(out,fieldnames=fields,delimiter='\t');w.writeheader()
  for x in r:
   kid=x['kmer_id'];ch=x['chromosome'];s=int(x['start_0based']);e=int(x['end_0based']);assert kid in fams and ch in CHRS and e-s==16
   br=match(broad,ch,s,e,list(broad[ch]));fi=match(fine,ch,s,e,['HSat3','HSat2','alphaSat','ct'])
   b=int(br!='outside');x.update(target_family=fams[kid],inside_broad_CenSat=b,fine_censat_family=fi,outside_broad_CenSat=1-b,outside_fine_annotation=int(fi=='outside'));w.writerow(x)
   counts[(kid,ch,'total')]+=1;counts[(kid,ch,'broad_inside')]+=b;counts[(kid,ch,'broad_outside')]+=1-b;counts[(kid,ch,'fine_'+fi)]+=1
   bins[(kid,ch,s//100000,1-b)]+=1;n+=1
 chrows=[]
 for kid in ids:
  for ch in CHRS:
   z=dict(kmer_id=kid,target_family=fams[kid],chromosome=ch,total_exact_hits=counts[(kid,ch,'total')],inside_broad_CenSat=counts[(kid,ch,'broad_inside')],outside_broad_CenSat=counts[(kid,ch,'broad_outside')]);z.update({f'fine_{label}_hits':counts[(kid,ch,'fine_'+label)] for label in LABELS});assert z['total_exact_hits']==z['inside_broad_CenSat']+z['outside_broad_CenSat']==sum(z[f'fine_{label}_hits'] for label in LABELS);chrows.append(z)
 write(OUT/'chm13_censat_by_chromosome.tsv',chrows,list(chrows[0]))
 totals=[]
 for kid in ids:
  rs=[r for r in chrows if r['kmer_id']==kid];z=dict(kmer_id=kid,target_family=fams[kid],total_exact_hits=sum(r['total_exact_hits'] for r in rs),inside_broad_CenSat=sum(r['inside_broad_CenSat'] for r in rs),outside_broad_CenSat=sum(r['outside_broad_CenSat'] for r in rs),chromosomes_with_hits=sum(r['total_exact_hits']>0 for r in rs));z.update({f'fine_{label}_hits':sum(r[f'fine_{label}_hits'] for r in rs) for label in LABELS});totals.append(z)
 assert sum(r['total_exact_hits'] for r in totals)==n
 write(OUT/'chm13_censat_totals.tsv',totals,list(totals[0]))
 ranked={r['kmer_id']:r for r in read(ROOT/'step4/ranked/mismatch_evaluated_ranked.tsv')}
 candidate=[]
 for r in totals:
  x=ranked[r['kmer_id']];candidate.append(dict(kmer_id=r['kmer_id'],canonical_target_5to3=x['canonical_target_5to3'],target_family=r['target_family'],candidate_type=x['candidate_type'],category_robustness_tier=x['category_robustness_tier'],final_score=x['final_score'],CHM13_total_exact_hits=r['total_exact_hits'],CHM13_inside_broad_CenSat=r['inside_broad_CenSat'],CHM13_outside_broad_CenSat=r['outside_broad_CenSat'],CHM13_outside_fraction=r['outside_broad_CenSat']/r['total_exact_hits'] if r['total_exact_hits'] else '',CHM13_fine_target_family_hits=r['fine_'+r['target_family']+'_hits'],CHM13_fine_other_family_hits=r['fine_'+('HSat3' if r['target_family']=='HSat2' else 'HSat2')+'_hits'],CHM13_fine_alphaSat_hits=r['fine_alphaSat_hits'],CHM13_fine_ct_hits=r['fine_ct_hits'],CHM13_outside_fine_annotation_hits=r['fine_outside_hits'],chromosomes_with_hits=r['chromosomes_with_hits']))
 candidate.sort(key=lambda z:(-z['CHM13_outside_broad_CenSat'],-z['CHM13_outside_fraction'] if z['CHM13_outside_fraction']!='' else 0,z['kmer_id']))
 write(OUT/'chm13_candidate_censat_summary.tsv',candidate,list(candidate[0]))
 binrows=[dict(kmer_id=kid,chromosome=ch,bin_100kb=b,bin_start_0based=b*100000,outside_broad_CenSat=outside,exact_hits=count) for (kid,ch,b,outside),count in bins.items()]
 write(OUT/'chm13_100kb_bins.tsv',sorted(binrows,key=lambda r:(ids.index(r['kmer_id']),CHRS.index(r['chromosome']),r['bin_100kb'],r['outside_broad_CenSat'])),list(binrows[0]))
 groups={r['hap']:r for r in read(ROOT/'step3/hap_groups.tsv')};haprows=[];group_total=[];group_chr=[]
 with np.load(ROOT/'step4/mismatch/combined_hap_counts.npz') as z:
  c=z['chromosome_counts'];den=z['chromosome_valid_starts'];av=z['chromosome_available'];haps=list(z['haplotypes']);codes=list(z['query_ids']);assert codes==ids and list(z['chromosomes'])==CHRS and len(haps)==574
  # Independently scanned fine HSat2/HSat3 coordinates must reproduce STEP4's CHM13 0-mm matrix per chromosome.
  for qi,kid in enumerate(ids):
   for ci,ch in enumerate(CHRS):
    for fi,fam in enumerate(['HSat2','HSat3']):assert counts[(kid,ch,'fine_'+fam)]==int(c[qi,0,ci,fi,0]),(kid,ch,fam)
  detail=OUT/'hap_chromosome_exact_counts.tsv.gz';fields2=['kmer_id','target_family','hap','group','group_status','chromosome','chromosome_available','target_family_exact_hits','target_family_valid_16mer_starts','target_family_density_per_M','HSat23_union_exact_hits','background_exact_hits','union_plus_background_exact_hits']
  with gzip.open(detail,'wt') as f:
   w=csv.DictWriter(f,fieldnames=fields2,delimiter='\t');w.writeheader()
   for qi,kid in enumerate(ids):
    fam=fams[kid];fi=['HSat2','HSat3'].index(fam);grpfield=fam+'_group';statusfield=fam+'_status';group_names=sorted({row[grpfield] for row in groups.values() if row[grpfield].startswith('G')})
    for hi,hap in enumerate(haps):
     group=groups[hap][grpfield] if hi else 'reference';status=groups[hap][statusfield];target=int(c[qi,hi,:,fi,0].sum());windows=int(den[hi,:,fi].sum());union=int(c[qi,hi,:,2,0].sum());bg=int(c[qi,hi,:,3,0].sum())
     haprows.append(dict(kmer_id=kid,target_family=fam,hap=hap,group=group,group_status=status,target_family_exact_hits=target,target_family_valid_16mer_starts=windows,target_family_density_per_M=target/windows*1e6 if windows else '',HSat23_union_exact_hits=union,background_exact_hits=bg,union_plus_background_exact_hits=union+bg))
     for ci,ch in enumerate(CHRS):
      n1=int(c[qi,hi,ci,fi,0]);d=int(den[hi,ci,fi]);u=int(c[qi,hi,ci,2,0]);b=int(c[qi,hi,ci,3,0]);w.writerow(dict(kmer_id=kid,target_family=fam,hap=hap,group=group,group_status=status,chromosome=ch,chromosome_available=int(av[hi,ci]),target_family_exact_hits=n1,target_family_valid_16mer_starts=d,target_family_density_per_M=n1/d*1e6 if d else '',HSat23_union_exact_hits=u,background_exact_hits=b,union_plus_background_exact_hits=u+b))
    for group in group_names+['NA']:
     hi=np.array([j for j,hap in enumerate(haps) if j>0 and groups[hap][grpfield]==group]);if_n=len(hi)
     assert if_n>0
     vals=c[qi,hi,:,fi,0].sum(axis=1).astype(float);d=den[hi,:,fi].sum(axis=1).astype(float);density=vals[d>0]/d[d>0]*1e6
     group_total.append(dict(kmer_id=kid,target_family=fam,group=group,hap_count=if_n,median_target_hits=float(np.median(vals)),mean_target_hits=float(vals.mean()),q25_target_hits=float(np.quantile(vals,.25)),q75_target_hits=float(np.quantile(vals,.75)),median_target_density_per_M=float(np.median(density)) if len(density) else '',measurable_haps=len(density)))
     for ci,ch in enumerate(CHRS):
      measurable=av[hi,ci]&(den[hi,ci,fi]>0);v=c[qi,hi[measurable],ci,fi,0].astype(float);dw=den[hi[measurable],ci,fi].astype(float)
      group_chr.append(dict(kmer_id=kid,target_family=fam,group=group,chromosome=ch,available_haps=int(av[hi,ci].sum()),measurable_haps=len(v),median_target_hits=float(np.median(v)) if len(v) else '',mean_target_hits=float(v.mean()) if len(v) else '',pooled_target_density_per_M=float(v.sum()/dw.sum()*1e6) if len(v) else '',presence_ge10_fraction=float((v>=10).mean()) if len(v) else ''))
 write(OUT/'hap_total_exact_counts.tsv.gz',haprows,list(haprows[0]));write(OUT/'group_total_summary.tsv',group_total,list(group_total[0]));write(OUT/'group_chromosome_summary.tsv',group_chr,list(group_chr[0]))
 summary={'status':'PASS','queries':72,'CHM13_exact_hits':n,'sequence_sets':574,'hap_chromosome_rows':72*574*24,'fine_HSat2_and_HSat3_intervals_match_STEP0':True,'all_72_CHM13_fine_family_counts_match_STEP4_matrix':True,'broad_CenSat_BED':str(BROAD),'broad_CenSat_BED_sha256':sha(BROAD),'fine_CenSat_BED':str(FINE),'fine_CenSat_BED_sha256':sha(FINE),'coordinate_system':'0-based half-open exact 16-mers, canonical and reverse complement','group_scope':'STEP3 regional common-core groups; NA means incomplete core'}
 (OUT/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))
if __name__=='__main__':main()
