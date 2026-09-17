#!/usr/bin/env python3
"""Overlay fine CenSat families and broad CenSat windows on saved CHM13 hits."""
import csv,gzip,json,hashlib,bisect
from collections import defaultdict
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'step4/top_kmer_landscape'
FINE=ROOT.parent/'VallePrep_v0.0.0/workdir/work1/bed_all/chm13v2.0.censat.bed'
BROAD=ROOT.parent/'VallePrep_v0.0.0/reference/chm13v2.0_censat_v2.1.bed'
KIDS=['k16_3e476af8','k16_a0d0562c','k16_4118ba0e']
CHRS=[f'chr{i}' for i in range(1,23)]+['chrX','chrY']
def hashfile(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def readbed(p):
 out=defaultdict(lambda:defaultdict(list))
 for line in p.read_text().splitlines():
  if not line or line.startswith('#'):continue
  x=line.split();out[x[0]][x[3]].append((int(x[1]),int(x[2])))
 return {ch:{label:(tuple(r[0] for r in sorted(rs)),tuple(r[1] for r in sorted(rs))) for label,rs in by.items()} for ch,by in out.items()}
def family(bed,ch,start,end):
 if ch not in bed:return 'outside'
 for label in ['HSat3','HSat2','alphaSat','ct']+list(bed[ch]):
  if label not in bed[ch]:continue
  starts,ends=bed[ch][label];j=bisect.bisect_right(starts,start)-1
  while j>=0:
   if end<=ends[j]:return label
   j-=1
 return 'outside'
def main():
 fine=readbed(FINE);broad=readbed(BROAD)
 assert set(fine)==set(CHRS) and set(broad)==set(CHRS)
 for bed in [fine,broad]:
  for ch,by in bed.items():
   for label,(starts,ends) in by.items():assert all(starts[i]<ends[i] and (i==0 or starts[i]>=starts[i-1]) for i in range(len(starts))),(ch,label,len(starts))
 hsat3={(ch,s,e) for ch,by in fine.items() for s,e in zip(*by.get('HSat3',((),())))}
 existing={(x[0],int(x[1]),int(x[2])) for line in (ROOT/'step0/regions/chm13v2.0/HSat3.bed').read_text().splitlines() if (x:=line.split())}
 assert hsat3==existing
 src=OUT/'chm13_all_exact_hits.tsv.gz';dest=OUT/'chm13_all_exact_hits_with_censat.tsv.gz';counts=defaultdict(int);fields=None
 with gzip.open(src,'rt') as inp,gzip.open(dest,'wt') as out:
  r=csv.DictReader(inp,delimiter='\t');fields=r.fieldnames+['broad_censat','fine_censat_family','outside_broad_censat','outside_fine_censat_annotation'];w=csv.DictWriter(out,fieldnames=fields,delimiter='\t');w.writeheader()
  for x in r:
   ch=x['chromosome'];s=int(x['start_0based']);e=int(x['end_0based']);b=family(broad,ch,s,e);f=family(fine,ch,s,e)
   x.update(broad_censat=int(b!='outside'),fine_censat_family=f,outside_broad_censat=int(b=='outside'),outside_fine_censat_annotation=int(f=='outside'))
   w.writerow(x);counts[(x['kmer_id'],ch,'total')]+=1;counts[(x['kmer_id'],ch,'broad_'+('inside' if b!='outside' else 'outside'))]+=1;counts[(x['kmer_id'],ch,'fine_'+f)]+=1
 rows=[]
 for kid in KIDS:
  for ch in CHRS:
   z=dict(kmer_id=kid,chromosome=ch,total_exact_hits=counts[(kid,ch,'total')],inside_broad_CenSat=counts[(kid,ch,'broad_inside')],outside_broad_CenSat=counts[(kid,ch,'broad_outside')])
   for f in ['HSat3','HSat2','alphaSat','ct','outside']:z['fine_'+f+'_hits']=counts[(kid,ch,'fine_'+f)]
   assert z['total_exact_hits']==z['inside_broad_CenSat']+z['outside_broad_CenSat']
   assert z['total_exact_hits']==sum(z['fine_'+f+'_hits'] for f in ['HSat3','HSat2','alphaSat','ct','outside'])
   rows.append(z)
 with (OUT/'chm13_censat_inside_outside_by_chromosome.tsv').open('w') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
 totals=[]
 for kid in KIDS:
  rs=[r for r in rows if r['kmer_id']==kid];totals.append(dict(kmer_id=kid,total_exact_hits=sum(r['total_exact_hits'] for r in rs),inside_broad_CenSat=sum(r['inside_broad_CenSat'] for r in rs),outside_broad_CenSat=sum(r['outside_broad_CenSat'] for r in rs),fine_HSat3_hits=sum(r['fine_HSat3_hits'] for r in rs),fine_HSat2_hits=sum(r['fine_HSat2_hits'] for r in rs),fine_alphaSat_hits=sum(r['fine_alphaSat_hits'] for r in rs),fine_ct_hits=sum(r['fine_ct_hits'] for r in rs),outside_fine_annotation_hits=sum(r['fine_outside_hits'] for r in rs)))
 with (OUT/'chm13_censat_inside_outside_totals.tsv').open('w') as f:
  w=csv.DictWriter(f,fieldnames=list(totals[0]),delimiter='\t');w.writeheader();w.writerows(totals)
 assert all(r['total_exact_hits']==r['inside_broad_CenSat']+r['outside_broad_CenSat'] for r in totals)
 assert all(r['total_exact_hits']==sum(r[k] for k in ['fine_HSat3_hits','fine_HSat2_hits','fine_alphaSat_hits','fine_ct_hits','outside_fine_annotation_hits']) for r in totals)
 summary={'status':'PASS','fine_censat_bed':str(FINE),'fine_censat_bed_sha256':hashfile(FINE),'broad_censat_bed':str(BROAD),'broad_censat_bed_sha256':hashfile(BROAD),'fine_HSat3_intervals_match_STEP0':len(hsat3),'coordinate_system':'0-based half-open; each exact 16-mer must be completely inside annotation interval','kmer_totals':totals}
 (OUT/'censat_overlay_summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))
if __name__=='__main__':main()
