#!/usr/bin/env python3
"""Jellyfish-based exact 16-mer discovery and candidate-restricted region counts."""
import argparse, csv, gzip, hashlib, json, math, os, subprocess, time
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
import numpy as np
from workflow_settings import settings
ROOT=Path(__file__).resolve().parents[1]; OUT=ROOT/'step1'; S0=ROOT/'step0'
JF=ROOT/'.local/bin/jellyfish'; PACK=ROOT/'scripts/step1_pack'
SAM=Path(settings()['samtools'])
REGIONS=['HSat2','HSat3','HSat23','background']
DT=np.dtype([('code','<u4'),('count','<u8')])
CONFIG={'k':16,'epsilon':1e-9,'min_log2_enrichment':2.0,'min_pooled_target_count':100,'min_hap_target_count':10,'reference_policy':'CHM13 excluded from pooled discovery, evaluated separately','candidate_rule':'pooled enrichment>=2 and target>=100 OR at least one nonreference hap enrichment>=2 and target>=10','region_order':REGIONS,'orientation_order':['canonical_forward','reverse_complement'],'software':'Jellyfish 2.3.1','samtools':str(SAM),'count_precision':'Jellyfish uint32 per hap; uint64 downstream', 'jellyfish_out_counter_bytes':4,'preliminary_rule':'nonreference pooled target>=100 OR nonreference maximum hap target>=10','background':'complement of existing HSat2/3 BED on eligible T2T chromosomes'}

def table(path):
 with open(path) as f:return list(csv.DictReader(f,delimiter='\t'))
MANIFEST=[r for r in table(S0/'manifest.tsv') if r['status']=='PASS']
STATS={}
for r in table(S0/'region_stats.tsv'):
 key=(r['hap'],r['region']);STATS[key]=STATS.get(key,0)+int(r['valid_16mer_starts'])

def fingerprint(r,region,seed=None):
 bed=S0/'regions'/r['hap']/(region+'.bed')
 st=Path(r['fasta']).stat()
 value=[r['fasta'],st.st_size,st.st_mtime_ns,hashlib.sha256(bed.read_bytes()).hexdigest(),CONFIG,seed]
 return hashlib.sha256(json.dumps(value,sort_keys=True).encode()).hexdigest()

def readbin(path):
 with gzip.open(path,'rb') as f:return np.frombuffer(f.read(),dtype=DT).copy()

def revcodes(x):
 x=x.astype(np.uint32,copy=True);out=np.zeros_like(x)
 for _ in range(16):out=(out<<np.uint32(2))|((x&np.uint32(3))^np.uint32(3));x>>=np.uint32(2)
 return out

def sequence(code):
 return ''.join('ACGT'[(int(code)>>(2*i))&3] for i in range(15,-1,-1))

def count_region(r,region,dest,canonical=True,seed=None,hashsize='8M'):
 hap=r['hap'];assert STATS[(hap,region)] < 2**32, 'Use safe counter backend for this input';work=OUT/'work'/hap;work.mkdir(parents=True,exist_ok=True)
 bed=S0/'regions'/hap/(region+'.bed');regions=work/(region+'.regions.txt')
 lines=[]
 for line in bed.read_text().splitlines():
  c,s,e,*_=line.split();lines.append(f'{c}:{int(s)+1}-{e}\n')
 regions.write_text(''.join(lines))
 if not lines:
  with gzip.open(dest,'wb') as f:pass
  return {'distinct':0,'count_sum':0,'seconds':0}
 database=work/(region+'.jf');log=OUT/'logs'/f'{hap}.{region}.log';timer=OUT/'logs'/f'{hap}.{region}.time.txt'
 counter=ROOT/'.local/bin/jellyfish-prefix' if seed and (ROOT/'.local/bin/jellyfish-prefix').exists() else JF
 cmd=[str(counter),'count','-m','16','-s',str(hashsize),'-t','2','-c','16','--out-counter-len','4','-L','1','-o',str(database)]
 if canonical:cmd.append('-C')
 if seed:cmd+=['--if',str(seed)]
 cmd.append('/dev/stdin')
 start=time.monotonic()
 with log.open('w') as err:
  producer=subprocess.Popen([str(SAM),'faidx','-n','1000000','-r',str(regions),r['fasta']],stdout=subprocess.PIPE,stderr=err)
  try:
   env=os.environ.copy()
   if counter.name=='jellyfish-prefix':env['SAT_ASO_PREFIX_FILTER']='1'
   counted=subprocess.run(['/usr/bin/time','-v','-o',str(timer)]+cmd,stdin=producer.stdout,stderr=err,env=env)
  finally:producer.stdout.close()
  pret=producer.wait()
  if counted.returncode or pret:raise RuntimeError(f'count failed {hap} {region}: {log}')
  dump=subprocess.Popen([str(JF),'dump','-c',str(database)],stdout=subprocess.PIPE,stderr=err)
  packed=subprocess.run([str(PACK),'pack',str(dest)],stdin=dump.stdout,stdout=subprocess.PIPE,stderr=err,text=True)
  dump.stdout.close();dret=dump.wait()
  if packed.returncode or dret:raise RuntimeError(f'dump failed {hap} {region}')
 distinct,total=map(int,packed.stdout.split())
 database.unlink();regions.unlink()
 return {'distinct':distinct,'count_sum':total,'seconds':time.monotonic()-start,'counter_binary':str(counter),'exact_prefix_gate':counter.name=='jellyfish-prefix'}

def discovery(r):
 hap=r['hap'];dest=OUT/'discovery'/(hap+'.bin.gz');meta=dest.with_suffix('.json');sig=fingerprint(r,'HSat23')
 if dest.exists() and meta.exists() and json.loads(meta.read_text()).get('fingerprint')==sig:return hap,'cached'
 tmp=dest.with_suffix('.tmp.gz')
 qc=count_region(r,'HSat23',tmp)
 expected=STATS[(hap,'HSat23')]
 if qc['count_sum']!=expected:raise RuntimeError(f'count sum {qc} != step0 {expected} for {hap}')
 tmp.replace(dest);qc.update(fingerprint=sig,expected_valid_starts=expected)
 meta.write_text(json.dumps(qc,indent=2)+'\n');return hap,qc

def make_pool():
 paths=[OUT/'discovery'/(r['hap']+'.bin.gz') for r in MANIFEST if r['sample']!='CHM13']
 if not all(p.exists() for p in paths):raise RuntimeError('discovery incomplete')
 listing=OUT/'discovery'/'nonreference_files.txt';listing.write_text(''.join(str(p)+'\n' for p in paths))
 res=subprocess.run([str(PACK),'pool',str(listing),str(OUT/'preliminary_candidates.u32')],capture_output=True,text=True,check=True)
 (OUT/'discovery_summary.json').write_text(res.stdout)
 codes=np.fromfile(OUT/'preliminary_candidates.u32',dtype='<u4');rc=revcodes(codes)
 with (OUT/'preliminary_candidates.both_strands.fa').open('w') as f:
  for i,(code,rev) in enumerate(zip(codes,rc)):
   f.write(f'>{i}_F\n{sequence(code)}\n')
   if code!=rev:f.write(f'>{i}_R\n{sequence(rev)}\n')
 bits=np.zeros(1<<18,dtype='<u8');prefixes=np.r_[codes,rc].astype(np.uint64)>>8
 np.bitwise_or.at(bits,prefixes>>6,np.left_shift(np.uint64(1),prefixes&63))
 bits.tofile(str(OUT/'preliminary_candidates.both_strands.fa')+'.prefix12.bitmap')
 print('POOL',res.stdout.strip(),flush=True)

def candidate_counts(r,codes,seed,seedhash):
 hap=r['hap'];dest=OUT/'counts'/(hap+'.npz');meta=OUT/'counts'/(hap+'.json')
 sig=hashlib.sha256(''.join(fingerprint(r,f,seedhash) for f in REGIONS).encode()).hexdigest()
 if dest.exists() and meta.exists() and json.loads(meta.read_text()).get('fingerprint')==sig:return hap,'cached'
 counts=np.zeros((len(codes),4,2),dtype=np.uint64);checks={}
 for ri,region in enumerate(REGIONS):
  rawfile=OUT/'work'/hap/(region+'.bin.gz');rawfile.parent.mkdir(parents=True,exist_ok=True)
  qc=count_region(r,region,rawfile,canonical=False,seed=seed if region=='background' else None,hashsize=max(8*1024**2,len(codes)*3) if region=='background' else '8M')
  raw=readbin(rawfile);rev=revcodes(raw['code']);can=np.minimum(raw['code'],rev);idx=np.searchsorted(codes,can)
  selected=idx<len(codes)
  selected[selected] &= codes[idx[selected]]==can[selected]
  if region=='background' and not selected.all():raise RuntimeError('unseeded background k-mer')
  orientation=(raw['code']>rev).astype(np.int8)
  counts[idx[selected],ri,orientation[selected]]=raw['count'][selected]
  if region!='background' and qc['count_sum']!=STATS[(hap,region)]:raise RuntimeError(f'full target count mismatch {hap} {region}')
  if qc['count_sum']>STATS[(hap,region)]:raise RuntimeError('more counts than valid windows')
  qc['expected_valid_starts']=STATS[(hap,region)];qc['candidate_count_sum']=int(counts[:,ri,:].sum());checks[region]=qc
  rawfile.unlink()
 # Independent canonical discovery must equal the sum of the two orientation counts.
 raw=readbin(OUT/'discovery'/(hap+'.bin.gz'));idx=np.searchsorted(codes,raw['code']);sel=idx<len(codes);sel[sel]&=codes[idx[sel]]==raw['code'][sel]
 expected=np.zeros(len(codes),dtype=np.uint64);expected[idx[sel]]=raw['count'][sel]
 if not np.array_equal(expected,counts[:,2,:].sum(axis=1)):raise RuntimeError('canonical/orientation mismatch '+hap)
 tmp=dest.with_suffix('.tmp.npz');np.savez_compressed(tmp,counts=counts);tmp.replace(dest)
 meta.write_text(json.dumps({'fingerprint':sig,'regions':checks,'canonical_orientation_agreement':True},indent=2)+'\n')
 return hap,{f:round(checks[f]['seconds'],1) for f in REGIONS}

def run_counts(jobs,workers):
 codes=np.fromfile(OUT/'preliminary_candidates.u32',dtype='<u4');seed=OUT/'preliminary_candidates.both_strands.fa'
 seedhash=hashlib.sha256((OUT/'preliminary_candidates.u32').read_bytes()).hexdigest()
 with ThreadPoolExecutor(max_workers=workers) as pool:
  for i,f in enumerate(as_completed([pool.submit(candidate_counts,r,codes,seed,seedhash) for r in jobs]),1):print(i,len(jobs),*f.result(),flush=True)

def enrichment(target,background,target_den,background_den):
 if target_den<=0 or background_den<=0:
  nan=np.full(len(target),np.nan);return nan,nan.copy(),nan.copy()
 td=target.astype(np.float64)/target_den;bd=background.astype(np.float64)/background_den
 return td,bd,np.log2((td+CONFIG['epsilon'])/(bd+CONFIG['epsilon']))

def score():
 codes=np.fromfile(OUT/'preliminary_candidates.u32',dtype='<u4');n=len(codes)
 pooled=np.zeros((n,4,2),dtype=np.uint64);den=np.zeros(4,dtype=np.uint64)
 prevalence=np.zeros(n,dtype=np.uint16);qualified=np.zeros(n,dtype=np.uint16)
 esum=np.zeros(n);emin=np.full(n,np.inf);emax=np.full(n,-np.inf);nscore=0;reference=None
 qcrows=[];haprows=[]
 for hi,r in enumerate(MANIFEST,1):
  hap=r['hap'];file=OUT/'counts'/(hap+'.npz')
  if not file.exists():raise RuntimeError('missing candidate counts '+hap)
  with np.load(file) as z:counts=z['counts']
  if counts.shape!=(n,4,2):raise RuntimeError('count shape mismatch')
  denominators=np.array([STATS[(hap,f)] for f in REGIONS],dtype=np.uint64)
  t=counts[:,2,:].sum(axis=1);b=counts[:,3,:].sum(axis=1)
  td,bd,e=enrichment(t,b,denominators[2],denominators[3])
  np.savez_compressed(OUT/'enrichment_by_hap'/(hap+'.npz'),target_density=td,background_density=bd,log2_enrichment=e,valid_starts=denominators)
  haprows.append([hap,r['sample']=='CHM13',int(denominators[2]),int(denominators[3]),int((t>0).sum()),int(((t>=10)&(e>=2)).sum())])
  if r['sample']=='CHM13':reference=(counts.copy(),denominators.copy(),e.copy())
  else:
   pooled+=counts;den+=denominators;prevalence+=(t>0).astype(np.uint16)
   qualified+=((t>=CONFIG['min_hap_target_count']) & (e>=CONFIG['min_log2_enrichment'])).astype(np.uint16)
   if denominators[2]>0 and denominators[3]>0:
    esum+=e;emin=np.minimum(emin,e);emax=np.maximum(emax,e);nscore+=1
  qc=json.loads((OUT/'counts'/(hap+'.json')).read_text())
  for f in REGIONS:
   v=qc['regions'][f];qcrows.append([hap,f,v['expected_valid_starts'],v['count_sum'],v['candidate_count_sum'],v['distinct'],v['seconds'],qc['canonical_orientation_agreement']])
  print('score',hi,len(MANIFEST),hap,flush=True)
 t=pooled[:,2,:].sum(axis=1);b=pooled[:,3,:].sum(axis=1)
 td,bd,e=enrichment(t,b,den[2],den[3])
 pooled_pass=(t>=CONFIG['min_pooled_target_count']) & (e>=CONFIG['min_log2_enrichment']);hap_pass=qualified>0;selected=pooled_pass|hap_pass
 np.savez_compressed(OUT/'pooled_counts_and_scores.npz',codes=codes,counts=pooled,valid_starts=den,target_density=td,background_density=bd,log2_enrichment=e,prevalence_haps=prevalence,qualifying_haps=qualified,mean_hap_log2_enrichment=esum/max(1,nscore),min_hap_log2_enrichment=emin,max_hap_log2_enrichment=emax,selected=selected,pooled_pass=pooled_pass,hap_pass=hap_pass)
 codes[selected].tofile(OUT/'candidate_codes.u32')
 header=['kmer_id','canonical_kmer','reverse_complement','pooled_target_count','pooled_background_count','pooled_target_density','pooled_background_density','pooled_log2_enrichment','target_prevalence_haps','target_prevalence_fraction','qualifying_haps','mean_hap_log2_enrichment','min_hap_log2_enrichment','max_hap_log2_enrichment','CHM13_target_count','CHM13_background_count','CHM13_log2_enrichment','selected','selection_reason']
 catalog_header=['kmer_id','canonical_kmer','forward_sequence','reverse_complement','orientation','reverse_complement_palindrome','GC_percent','max_homopolymer','base_entropy_bits']
 with gzip.open(OUT/'enrichment_summary.tsv.gz','wt') as summary, open(OUT/'candidate_kmers.tsv','w') as candidates, gzip.open(OUT/'kmer_catalog.tsv.gz','wt') as catalog:
  sw=csv.writer(summary,delimiter='\t');cw=csv.writer(candidates,delimiter='\t');kw=csv.writer(catalog,delimiter='\t')
  sw.writerow(header);cw.writerow(header+catalog_header[4:]);kw.writerow(catalog_header)
  for i,code in enumerate(codes):
   seq=sequence(code);reverse=seq.translate(str.maketrans('ACGT','TGCA'))[::-1];kid=f'k16_{int(code):08x}'
   run=maxrun=1
   for j in range(1,16):run=run+1 if seq[j]==seq[j-1] else 1;maxrun=max(maxrun,run)
   freqs=[seq.count(base)/16 for base in 'ACGT'];entropy=-sum(p*math.log2(p) for p in freqs if p)
   features=['canonical',seq==reverse,100*(seq.count('G')+seq.count('C'))/16,maxrun,entropy]
   kw.writerow([kid,seq,seq,reverse]+features)
   reason='pooled_and_hap' if pooled_pass[i] and hap_pass[i] else ('pooled' if pooled_pass[i] else ('hap' if hap_pass[i] else 'not_selected'))
   refvals=[int(reference[0][i,2,:].sum()),int(reference[0][i,3,:].sum()),reference[2][i]] if reference else ['NA']*3
   row=[kid,seq,reverse,int(t[i]),int(b[i]),td[i],bd[i],e[i],int(prevalence[i]),float(prevalence[i])/max(1,nscore),int(qualified[i]),esum[i]/max(1,nscore),emin[i],emax[i]]+refvals+[bool(selected[i]),reason]
   sw.writerow(row)
   if selected[i]:cw.writerow(row+features)
 with (OUT/'qc_summary.tsv').open('w') as f:
  w=csv.writer(f,delimiter='\t');w.writerow(['hap','region','valid_16mer_starts','all_or_selected_count_sum','preliminary_candidate_count_sum','distinct_counted_kmers','seconds','canonical_orientation_agreement']);w.writerows(qcrows)
 with (OUT/'hap_summary.tsv').open('w') as f:
  w=csv.writer(f,delimiter='\t');w.writerow(['hap','is_reference','target_valid_starts','background_valid_starts','preliminary_kmers_present_in_target','locally_qualifying_kmers']);w.writerows(haprows)
 thresholds=[]
 for threshold in [1,2,3,4,5]:
  for mincount in [10,100,1000]:thresholds.append([threshold,mincount,int(((e>=threshold)&(t>=mincount)).sum())])
 with (OUT/'pooled_threshold_sensitivity.tsv').open('w') as f:
  w=csv.writer(f,delimiter='\t');w.writerow(['min_log2_enrichment','min_pooled_target_count','pooled_candidates']);w.writerows(thresholds)
 summary={'sequence_sets':len(MANIFEST),'nonreference_haps':nscore,'preliminary_candidates':n,'selected_candidates':int(selected.sum()),'pooled_pass':int(pooled_pass.sum()),'hap_pass':int(hap_pass.sum()),'hap_only':int((hap_pass&~pooled_pass).sum()),'pooled_only':int((pooled_pass&~hap_pass).sum()),'both':int((pooled_pass&hap_pass).sum()),'pooled_target_valid_starts':int(den[2]),'pooled_background_valid_starts':int(den[3]),'epsilon':CONFIG['epsilon']}
 (OUT/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
 print('SUMMARY',summary,flush=True)


def main():
 ap=argparse.ArgumentParser();ap.add_argument('--phase',choices=['discovery','counts','score'],required=True);ap.add_argument('--workers',type=int,default=8);ap.add_argument('--pilot',type=int,default=0);args=ap.parse_args()
 OUT.mkdir(exist_ok=True)
 for sub in ['discovery','counts','logs','plots','enrichment_by_hap','work']:(OUT/sub).mkdir(exist_ok=True)
 (OUT/'config.json').write_text(json.dumps(CONFIG,indent=2)+'\n')
 jobs=MANIFEST[:args.pilot] if args.pilot else MANIFEST
 if args.phase=='discovery':
  with ThreadPoolExecutor(max_workers=args.workers) as pool:
   for i,future in enumerate(as_completed([pool.submit(discovery,r) for r in jobs]),1):print(i,len(jobs),*future.result(),flush=True)
  if not args.pilot:make_pool()
 elif args.phase=='counts':run_counts(jobs,args.workers)
 else:score()
if __name__=='__main__':main()
