#!/usr/bin/env python3
"""Chromosome-resolved exact candidate counts, with STEP1 reconciliation."""
import argparse,csv,hashlib,json,os,subprocess,time
from concurrent.futures import ThreadPoolExecutor,as_completed
from pathlib import Path
import numpy as np
from step1_pipeline import ROOT,S0,SAM,MANIFEST,sequence
from workflow_settings import settings

OUT=ROOT/'step2';S1=ROOT/'step1';REGIONS=['HSat2','HSat3','HSat23','ambiguous']
SELECT=settings()['selection']
SCANNER=ROOT/'scripts/step2_count';SOURCE=S1/'candidate_kmers.tsv'
ROWS=list(csv.DictReader((S0/'region_stats.tsv').open(),delimiter='\t'))
STATS={(r['hap'],r['chromosome'],r['region']):r for r in ROWS}
def chromkey(c):
 s=c.removeprefix('chr');return (int(s) if s.isdigit() else {'X':23,'Y':24}.get(s,99),c)
CHROMS=sorted({r['chromosome'] for r in ROWS},key=chromkey)
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def write_tsv(path,header,rows):
 with path.open('w') as f:
  w=csv.writer(f,delimiter='\t');w.writerow(header);w.writerows(rows)
def setup():
 for sub in ['counts','work','logs','matrices','plots']:(OUT/sub).mkdir(parents=True,exist_ok=True)
 src=list(csv.DictReader(SOURCE.open(),delimiter='\t'))
 codes=np.array([int(r['kmer_id'].split('_')[1],16) for r in src],dtype='<u4')
 assert len(codes)>0 and np.all(codes[1:]>codes[:-1])
 assert all(r['selected']=='True' and r['selection_reason'] in ('pooled','hap','pooled_and_hap') for r in src)
 candidate_file=OUT/'candidate_codes.u32'
 if candidate_file.exists():assert candidate_file.read_bytes()==codes.tobytes()
 else:codes.tofile(candidate_file)
 broad_min=int(np.ceil(SELECT['broad_hap_fraction']*sum(r['sample']!='CHM13' for r in MANIFEST)))
 write_tsv(OUT/'candidates.tsv',['row_index','kmer_id','canonical_kmer','reverse_complement','pooled_selected_threshold','pooled_and_broad_threshold'],((i,r['kmer_id'],r['canonical_kmer'],r['reverse_complement'],int(r['selection_reason'] in ('pooled','pooled_and_hap')),int(r['selection_reason'] in ('pooled','pooled_and_hap') and int(r['qualifying_haps'])>=broad_min)) for i,r in enumerate(src)))
 write_tsv(OUT/'haplotypes.tsv',['column_index','hap','is_reference','fasta'],((i,r['hap'],int(r['sample']=='CHM13'),r['fasta']) for i,r in enumerate(MANIFEST)))
 config={'k':16,'candidate_source':str(SOURCE),'selected_E':SELECT['selected_E'],'selection':SELECT,'candidate_sha256':sha(SOURCE),'candidate_count':len(codes),'sequence_sets':len(MANIFEST),'nonreference_haps':sum(r['sample']!='CHM13' for r in MANIFEST),'regions':REGIONS,'chromosomes':CHROMS,'strand_order':['canonical_forward','reverse_complement'],'missing_chromosomes':'absent from shard; availability mask false, not biological zero','counter_sha256':sha(SCANNER),'step0_complete_sha256':sha(S0/'COMPLETE.json'),'step1_complete_sha256':sha(S1/'COMPLETE.json')}
 (OUT/'config.json').write_text(json.dumps(config,indent=2)+'\n')
 return codes,config
def count_hap(r,codes,config):
 h=r['hap'];chrs=[c for c in CHROMS if (h,c,'HSat23') in STATS];nc=len(chrs)
 st=Path(r['fasta']).stat()
 sig=hashlib.sha256(json.dumps([config,st.st_size,st.st_mtime_ns,{f:sha(S0/'regions'/h/(f+'.bed')) for f in REGIONS},sha(S1/'counts'/(h+'.json'))],sort_keys=True).encode()).hexdigest()
 dest=OUT/'counts'/(h+'.npz');meta=dest.with_suffix('.json')
 if dest.exists() and meta.exists() and json.loads(meta.read_text())['fingerprint']==sig:return h,'cached'
 work=OUT/'work'/h;work.mkdir(exist_ok=True);(work/'chromosomes.txt').write_text(''.join(c+'\n' for c in chrs))
 result=np.zeros((len(codes),nc,4,2),dtype=np.uint64);den=np.zeros((nc,4),dtype=np.uint64);bp=den.copy();qc={};started=time.monotonic()
 for ri,region in enumerate(REGIONS):
  bed=S0/'regions'/h/(region+'.bed');lines=[]
  for line in bed.read_text().splitlines():
   c,s,e,*_=line.split();lines.append(f'{c}:{int(s)+1}-{e}\n')
  req=work/'regions.txt';req.write_text(''.join(lines));raw=work/'counts.u64';stats=work/'stats.tsv'
  command=[str(SCANNER),str(OUT/'candidate_codes.u32'),str(work/'chromosomes.txt'),str(raw),str(stats)]
  with (OUT/'logs'/f'{h}.{region}.log').open('w') as err:
   if lines:
    producer=subprocess.Popen([str(SAM),'faidx','-n','1000000','-r',str(req),r['fasta']],stdout=subprocess.PIPE,stderr=err)
    try: counted=subprocess.run(command,stdin=producer.stdout,stderr=err)
    finally:producer.stdout.close()
    pret=producer.wait()
    if pret or counted.returncode:raise RuntimeError(f'count failure {h} {region}')
   else:subprocess.run(command,input=b'',stderr=err,check=True)
  array=np.fromfile(raw,dtype='<u8').reshape(len(codes),nc,2);result[:,:,ri,:]=array
  scan=list(csv.DictReader(stats.open(),delimiter='\t'))
  for ci,row in enumerate(scan):
   assert row['chromosome']==chrs[ci]
   expected=STATS[(h,chrs[ci],region)]
   for key in ['bp','valid_16mer_starts']:assert int(row[key])==int(expected[key]),(h,region,chrs[ci],key)
   assert int(row['records'])==int(expected['interval_count'])
   den[ci,ri]=int(row['valid_16mer_starts']);bp[ci,ri]=int(row['bp'])
   assert int(array[:,ci,:].sum())<=int(den[ci,ri])
  qc[region]={'valid_windows':int(den[:,ri].sum()),'selected_count':int(array.sum())}
  raw.unlink()
 # Compare both orientations independently, for every candidate and target region.
 allcodes=np.fromfile(S1/'preliminary_candidates.u32',dtype='<u4');idx=np.searchsorted(allcodes,codes);assert np.array_equal(allcodes[idx],codes)
 with np.load(S1/'counts'/(h+'.npz')) as z:previous=z['counts'][idx,:3,:]
 assert np.array_equal(result[:,:,:3,:].sum(axis=1),previous),f'STEP1 mismatch: {h}'
 # Family windows are a subset of union windows; residual includes junction-spanning windows.
 assert np.all(result[:,:,2,:]>=result[:,:,0,:]+result[:,:,1,:]+result[:,:,3,:])
 tmp=dest.with_suffix('.tmp.npz');np.savez_compressed(tmp,counts=result,chromosomes=np.array(chrs),regions=np.array(REGIONS),valid_starts=den,bp=bp);tmp.replace(dest)
 meta.write_text(json.dumps({'fingerprint':sig,'hap':h,'chromosomes':chrs,'regions':qc,'step1_both_orientations_match':True,'step0_per_chromosome_windows_match':True,'seconds':time.monotonic()-started},indent=2)+'\n')
 return h,round(time.monotonic()-started,1)
def run(workers,pilot):
 codes,config=setup();jobs=MANIFEST[:pilot] if pilot else MANIFEST
 with ThreadPoolExecutor(max_workers=workers) as pool:
  futures=[pool.submit(count_hap,r,codes,config) for r in jobs]
  for i,f in enumerate(as_completed(futures),1):print('count',i,len(jobs),*f.result(),flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--workers',type=int,default=6);p.add_argument('--pilot',type=int,default=0);a=p.parse_args();run(a.workers,a.pilot)
