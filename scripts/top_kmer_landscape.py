#!/usr/bin/env python3
"""CHM13 genomic coordinates and group/chromosome exact counts for top HSat3 targets."""
import csv,gzip,json,bisect
from collections import defaultdict
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'step4/top_kmer_landscape';OUT.mkdir(parents=True,exist_ok=True)
KIDS=[]
def read(path):
 with open(path) as f:return list(csv.DictReader(f,delimiter='\t'))
def write(path,rows,fields):
 with (gzip.open(path,'wt') if str(path).endswith('.gz') else open(path,'w')) as f:
  w=csv.DictWriter(f,fieldnames=fields,delimiter='\t');w.writeheader();w.writerows(rows)
def rc(s):return s.translate(str.maketrans('ACGT','TGCA'))[::-1]
def fasta(path):
 with gzip.open(path,'rb') as f:
  name=None;parts=[]
  for line in f:
   if line.startswith(b'>'):
    if name is not None:yield name,b''.join(parts).upper()
    name=line[1:].split()[0].decode();parts=[]
   else:parts.append(line.strip())
  if name is not None:yield name,b''.join(parts).upper()
def intervals(path):
 by=defaultdict(list)
 for line in path.read_text().splitlines():
  if not line:continue
  ch,s,e,*_=line.split();by[ch].append((int(s),int(e)))
 return {ch:(np.array([x[0] for x in xs],dtype=np.int64),np.array([x[1] for x in xs],dtype=np.int64)) for ch,xs in by.items()}
def contains(reg,ch,p):
 if ch not in reg:return False
 starts,ends=reg[ch];j=bisect.bisect_right(starts,p)-1
 return j>=0 and p+16<=ends[j]
def main():
 global KIDS
 ranks=read(ROOT/'step4/ranked/mismatch_evaluated_ranked.tsv')
 KIDS=[r['kmer_id'] for r in ranks if r['target_family']=='HSat3' and r['candidate_type']=='pan'][:3]
 if len(KIDS)!=3:raise RuntimeError('fewer than three evaluated pan-HSat3 candidates')
 write(OUT/'selected_kmers.tsv',([{'kmer_id':kid} for kid in KIDS]),['kmer_id'])
 selected={r['kmer_id']:r for r in ranks if r['kmer_id'] in KIDS};assert len(selected)==3
 groups={(r['hap'],r['chromosome']):r for r in read(ROOT/'step3/hap_chromosome_groups.tsv') if r['family']=='HSat3'}
 manifest=read(ROOT/'step0/manifest.tsv');ref=next(r for r in manifest if r['sample']=='CHM13');fai={r.split()[0]:int(r.split()[1]) for r in Path(ref['fai']).read_text().splitlines()}
 chrs=[f'chr{i}' for i in range(1,23)]+['chrX','chrY']; assert set(chrs)==set(fai)
 reg={f:intervals(ROOT/'step0/regions/chm13v2.0'/f'{f}.bed') for f in ['HSat2','HSat3','HSat23']}
 motifs={kid:(selected[kid]['canonical_target_5to3'].encode(),selected[kid]['reverse_complement_target_5to3'].encode()) for kid in KIDS}
 hitfields=['kmer_id','canonical_kmer','chromosome','start_0based','end_0based','strand_of_canonical_target','region','bin_1Mb']; hitcounts=defaultdict(int);bins=defaultdict(int)
 with gzip.open(OUT/'chm13_all_exact_hits.tsv.gz','wt') as fh:
  w=csv.DictWriter(fh,fieldnames=hitfields,delimiter='\t');w.writeheader();seen=set()
  for ch,seq in fasta(ref['fasta']):
   if ch not in fai:continue
   assert len(seq)==fai[ch];seen.add(ch)
   for kid,(forward,reverse) in motifs.items():
    for strand,motif in [('+',forward),('-',reverse)]:
     start=0
     while True:
      p=seq.find(motif,start)
      if p<0:break
      region='HSat3' if contains(reg['HSat3'],ch,p) else ('HSat2' if contains(reg['HSat2'],ch,p) else ('HSat23_ambiguous' if contains(reg['HSat23'],ch,p) else 'background_or_boundary'))
      w.writerow(dict(kmer_id=kid,canonical_kmer=selected[kid]['canonical_target_5to3'],chromosome=ch,start_0based=p,end_0based=p+16,strand_of_canonical_target=strand,region=region,bin_1Mb=p//1000000))
      hitcounts[(kid,ch,region,strand)]+=1;bins[(kid,ch,p//1000000,region)]+=1;start=p+1
   print('scanned',ch,len(seq),flush=True)
 assert seen==set(chrs)
 chrows=[]
 for kid in KIDS:
  for ch in chrs:
   row=dict(kmer_id=kid,chromosome=ch,chromosome_length_bp=fai[ch])
   for region in ['HSat3','HSat2','HSat23_ambiguous','background_or_boundary']:
    row[region+'_hits']=sum(hitcounts[(kid,ch,region,strand)] for strand in ['+','-'])
   row['total_hits']=sum(row[f+'_hits'] for f in ['HSat3','HSat2','HSat23_ambiguous','background_or_boundary'])
   row['forward_hits']=sum(hitcounts[(kid,ch,r,'+')] for r in ['HSat3','HSat2','HSat23_ambiguous','background_or_boundary'])
   row['reverse_hits']=sum(hitcounts[(kid,ch,r,'-')] for r in ['HSat3','HSat2','HSat23_ambiguous','background_or_boundary'])
   chrows.append(row)
 write(OUT/'chm13_chromosome_summary.tsv',chrows,list(chrows[0]))
 binrows=[dict(kmer_id=kid,chromosome=ch,bin_1Mb=b,bin_start_0based=b*1000000,bin_end_0based=min((b+1)*1000000,fai[ch]),region=r,exact_hits=n) for (kid,ch,b,r),n in bins.items()]
 write(OUT/'chm13_1Mb_bins.tsv',sorted(binrows,key=lambda x:(KIDS.index(x['kmer_id']),chrs.index(x['chromosome']),x['bin_1Mb'],x['region'])),list(binrows[0]))

 with np.load(ROOT/'step4/mismatch/combined_hap_counts.npz') as z:
  counts=z['chromosome_counts'];den=z['chromosome_valid_starts'];av=z['chromosome_available'];ids=list(z['query_ids']);haps=list(z['haplotypes']);chroms=list(z['chromosomes'])
  assert chroms==chrs and haps[0]=='chm13v2.0';assert len(haps)==574
  rows=[];haptot=[];groupchrom=[]
  for kid in KIDS:
   qi=ids.index(kid)
   for hi,hap in enumerate(haps):
    t=0;d=0;allregions=0
    for ci,ch in enumerate(chrs):
     gr=groups[(hap,ch)];n=int(counts[qi,hi,ci,1,0]);windows=int(den[hi,ci,1]);union=int(counts[qi,hi,ci,2,0]);bg=int(counts[qi,hi,ci,3,0]);available=bool(av[hi,ci]);t+=n;d+=windows;allregions+=union+bg
     rows.append(dict(kmer_id=kid,hap=hap,HSat3_group=gr['group'],group_status=gr['status'],chromosome=ch,chromosome_available=int(available),HSat3_exact_hits=n,HSat3_valid_16mer_starts=windows,HSat3_density_per_M=n/windows*1e6 if windows else '',HSat23_union_exact_hits=union,background_exact_hits=bg,union_plus_background_exact_hits=union+bg))
    haptot.append(dict(kmer_id=kid,hap=hap,HSat3_exact_hits=t,HSat3_valid_16mer_starts=d,HSat3_density_per_M=t/d*1e6 if d else '',union_plus_background_exact_hits=allregions))
   for ci,ch in enumerate(chrs):
    names=sorted({groups[(hap,ch)]['group'] for hap in haps[1:] if groups[(hap,ch)]['group']!='NA'})
    for group in names:
     chosen=np.array([hi for hi,hap in enumerate(haps) if hi>0 and groups[(hap,ch)]['group']==group],dtype=int)
     v=counts[qi,chosen,ci,1,0].astype(float);dw=den[chosen,ci,1].astype(float)
     assert len(chosen)>0 and np.all(dw>0)
     groupchrom.append(dict(kmer_id=kid,HSat3_group=group,chromosome=ch,HSat3_measurable_haps=len(chosen),median_HSat3_hits=float(np.median(v)),mean_HSat3_hits=float(v.mean()),pooled_HSat3_density_per_M=float(v.sum()/dw.sum()*1e6),presence_ge10_fraction=float((v>=10).mean())))
 write(OUT/'hap_chromosome_exact_counts.tsv.gz',rows,list(rows[0]));write(OUT/'hap_total_exact_counts.tsv',haptot,list(haptot[0]))
 write(OUT/'group_chromosome_summary.tsv',groupchrom,list(groupchrom[0]) if groupchrom else ['kmer_id','HSat3_group','chromosome','HSat3_measurable_haps','median_HSat3_hits','mean_HSat3_hits','pooled_HSat3_density_per_M','presence_ge10_fraction'])
 # The independent CHM13 FASTA scan should agree for hits fully inside STEP0's HSat3 intervals.
 for kid in KIDS:
  for ch in chrs:
   direct=next(r['HSat3_hits'] for r in chrows if r['kmer_id']==kid and r['chromosome']==ch)
   matrix=next(r['HSat3_exact_hits'] for r in rows if r['kmer_id']==kid and r['hap']=='chm13v2.0' and r['chromosome']==ch)
   assert direct==matrix,(kid,ch,direct,matrix)
 summary={'status':'PASS','kmer_ids':KIDS,'chm13_filtered_chromosomes':len(chrs),'sequence_sets':len(haps),'nonreference_haps':len(haps)-1,'CHM13_HSat3_Fasta_hits_match_STEP4_matrix':True,'coordinate_system':'0-based half-open, exact 16-mer both orientations','group_scope':'HSat3 independent chromosome groups; NA means missing or unsupported'}
 (OUT/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))
if __name__=='__main__':main()
