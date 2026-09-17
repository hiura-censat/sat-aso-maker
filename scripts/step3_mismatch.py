#!/usr/bin/env python3
import argparse,csv,hashlib,json,subprocess,time
from concurrent.futures import ThreadPoolExecutor,as_completed
import numpy as np
from step2_pipeline import S0,S1,SAM,MANIFEST,STATS,CHROMS,sha
from step3_prepare import ROOT,OUT,S2
M=OUT/'mismatch';REGIONS=['HSat2','HSat3','HSat23','background'];SCANNER=ROOT/'scripts/step3_count'
def setup():
 for sub in ['counts','variants','work','logs']:(M/sub).mkdir(exist_ok=True)
 with np.load(M/'variant_mapping.npz') as z:records=z['records'];qcodes=z['query_codes'];codes=z['variant_codes']
 positions=[]
 for vi,qi,d,m1,m2 in records:
  for pos in range(16):
   weight=((int(m1)>>pos)&1)/2+((int(m2)>>pos)&1)/2
   if weight:positions.append((int(qi)*48+int(d)*16+pos,int(vi),weight))
 p=np.array(positions);signature=sha(M/'variant_mapping.npz')+sha(SCANNER)
 return records,qcodes,codes,p,signature
def count_hap(r,shared):
 records,qcodes,codes,pos,sigbase=shared;h=r['hap'];chrs=[c for c in CHROMS if (h,c,'HSat23') in STATS];nc=len(chrs)
 st=__import__('pathlib').Path(r['fasta']).stat();sig=hashlib.sha256(json.dumps([sigbase,st.st_size,st.st_mtime_ns,[sha(S0/'regions'/h/(f+'.bed')) for f in REGIONS]]).encode()).hexdigest()
 dest=M/'counts'/(h+'.npz');meta=dest.with_suffix('.json')
 if dest.exists() and meta.exists() and json.loads(meta.read_text())['fingerprint']==sig:return h,'cached'
 work=M/'work'/h;work.mkdir(exist_ok=True);(work/'chromosomes.txt').write_text(''.join(c+'\n' for c in chrs));started=time.monotonic()
 rawcounts=np.zeros((len(codes),nc,4,2),dtype=np.uint64);den=np.zeros((nc,4),dtype=np.uint64);timing={}
 for ri,region in enumerate(REGIONS):
  start=time.monotonic();lines=[]
  for line in (S0/'regions'/h/(region+'.bed')).read_text().splitlines():
   c,s,e,*_=line.split();lines.append(f'{c}:{int(s)+1}-{e}\n')
  req=work/'regions.txt';req.write_text(''.join(lines));binary=work/'counts.u64';stats=work/'stats.tsv'
  command=[str(SCANNER),str(M/'variant_codes.u32'),str(work/'chromosomes.txt'),str(binary),str(stats)]
  with (M/'logs'/f'{h}.{region}.log').open('w') as err:
   if lines:
    producer=subprocess.Popen([str(SAM),'faidx','-n','1000000','-r',str(req),r['fasta']],stdout=subprocess.PIPE,stderr=err)
    try:counted=subprocess.run(command,stdin=producer.stdout,stderr=err)
    finally:producer.stdout.close()
    pret=producer.wait()
    if pret or counted.returncode:raise RuntimeError(f'failed {h} {region}')
   else:subprocess.run(command,input=b'',stderr=err,check=True)
  arr=np.fromfile(binary,dtype='<u8').reshape(len(codes),nc,2);rawcounts[:,:,ri,:]=arr
  for ci,row in enumerate(csv.DictReader(stats.open(),delimiter='\t')):
   assert row['chromosome']==chrs[ci];old=STATS[(h,chrs[ci],region)]
   for key in ['bp','valid_16mer_starts']:assert int(row[key])==int(old[key]),(h,region,chrs[ci],key)
   assert int(row['records'])==int(old['interval_count']);den[ci,ri]=int(row['valid_16mer_starts']);assert int(arr[:,ci,:].sum())<=den[ci,ri]
  binary.unlink();timing[region]=round(time.monotonic()-start,2)
 # Distance zero must reproduce the independently counted STEP1 result including background.
 allcodes=np.fromfile(S1/'preliminary_candidates.u32',dtype='<u4');oldidx=np.searchsorted(allcodes,qcodes);queryidx=np.searchsorted(codes,qcodes)
 assert np.array_equal(allcodes[oldidx],qcodes) and np.array_equal(codes[queryidx],qcodes)
 with np.load(S1/'counts'/(h+'.npz')) as z:old=z['counts'][oldidx]
 assert np.array_equal(rawcounts[queryidx].sum(axis=1),old),f'exact mismatch {h}'
 canonical=rawcounts.sum(axis=3);acc=np.zeros((len(qcodes)*3,nc,4),dtype=np.uint64);np.add.at(acc,records[:,1]*3+records[:,2],canonical[records[:,0]])
 result=acc.reshape(len(qcodes),3,nc,4).transpose(0,2,3,1);assert np.all(result.sum(axis=3)<=den[None,:,:])
 positions=np.zeros((len(qcodes)*48,4));np.add.at(positions,pos[:,0].astype(int),canonical.sum(axis=1)[pos[:,1].astype(int)]*pos[:,2,None]);positions=positions.reshape(len(qcodes),3,16,4).transpose(0,3,1,2)
 for d in range(3):assert np.allclose(positions[:,:,d,:].sum(axis=2),result[:,:,:,d].sum(axis=1)*d)
 tmp=dest.with_suffix('.tmp.npz');np.savez_compressed(tmp,counts=result,position_counts=positions,chromosomes=np.array(chrs),regions=np.array(REGIONS),valid_starts=den);tmp.replace(dest)
 np.savez_compressed(M/'variants'/(h+'.npz'),counts=canonical,chromosomes=np.array(chrs),regions=np.array(REGIONS))
 meta.write_text(json.dumps({'hap':h,'fingerprint':sig,'zero_mismatch_matches_step1':True,'all_valid_windows_match_step0':True,'mismatch_position_accounting':True,'seconds':round(time.monotonic()-started,2),'region_seconds':timing},indent=2)+'\n');return h,timing
def run(workers,pilot):
 shared=setup();jobs=MANIFEST[:pilot] if pilot else MANIFEST
 with ThreadPoolExecutor(max_workers=workers) as pool:
  for i,f in enumerate(as_completed([pool.submit(count_hap,r,shared) for r in jobs]),1):print('mismatch',i,len(jobs),*f.result(),flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--workers',type=int,default=8);p.add_argument('--pilot',type=int,default=0);a=p.parse_args();run(a.workers,a.pilot)
