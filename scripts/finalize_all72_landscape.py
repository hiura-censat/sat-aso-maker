#!/usr/bin/env python3
"""Audit all 72 landscape outputs and build a navigable plot index."""
import csv,gzip,hashlib,json,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'step4/all72_landscape'
def read(path):
 with (gzip.open(path,'rt') if str(path).endswith('.gz') else open(path)) as f:return list(csv.DictReader(f,delimiter='\t'))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def linecount(path):
 with gzip.open(path,'rt') as f:return sum(1 for _ in f)-1
def main():
 source=sha(ROOT/'step4/COMPLETE.json');q=read(ROOT/'step4/mismatch/queries.tsv');ids=[r['kmer_id'] for r in q];assert ids and len(set(ids))==len(ids)
 s=json.loads((OUT/'summary.json').read_text());assert s['status']=='PASS' and s['queries']==len(ids) and s['all_72_CHM13_fine_family_counts_match_STEP4_matrix']
 totals=read(OUT/'chm13_censat_totals.tsv');detail=read(OUT/'chm13_censat_by_chromosome.tsv');candidates=read(OUT/'chm13_candidate_censat_summary.tsv');groups=read(OUT/'group_total_summary.tsv');gc=read(OUT/'group_chromosome_summary.tsv')
 assert len(totals)==len(ids) and len(detail)==len(ids)*24 and len(candidates)==len(ids)
 assert {r['kmer_id'] for r in totals}==set(ids)
 assert all(int(x['total_exact_hits'])==int(x['inside_broad_CenSat'])+int(x['outside_broad_CenSat']) for x in totals)
 assert all(int(x['total_exact_hits'])==sum(int(x[f'fine_{r}_hits']) for r in ['HSat3','HSat2','alphaSat','ct','outside']) for x in totals)
 assert linecount(OUT/'chm13_all72_exact_hits_with_censat.tsv.gz')==s['CHM13_exact_hits']
 nh=s['sequence_sets'];assert linecount(OUT/'hap_chromosome_exact_counts.tsv.gz')==len(ids)*nh*24
 assert linecount(OUT/'hap_total_exact_counts.tsv.gz')==len(ids)*nh
 assert len(groups)>0 and len(gc)==len(groups)*24
 # Independent old single-query Python scan must agree with the new one-pass C++ scanner for the original top three.
 old=ROOT/'step4/top_kmer_landscape/chm13_all_exact_hits.tsv.gz';original={}
 with gzip.open(old,'rt') as f:
  for x in csv.DictReader(f,delimiter='\t'):original.setdefault(x['kmer_id'],set()).add((x['chromosome'],x['start_0based'],x['strand_of_canonical_target']))
 recovered={kid:set() for kid in original}
 with gzip.open(OUT/'chm13_all72_exact_hits_with_censat.tsv.gz','rt') as f:
  for x in csv.DictReader(f,delimiter='\t'):
   if x['kmer_id'] in recovered:recovered[x['kmer_id']].add((x['chromosome'],x['start_0based'],x['strand_of_canonical_target']))
 assert recovered==original and len(original)==3
 plots=sorted((OUT/'plots').glob('*_landscape.png'));assert len(plots)==len(ids) and {p.name.removesuffix('_landscape.png') for p in plots}==set(ids)
 assert all(p.stat().st_size>100000 for p in plots)
 zooms=sorted((OUT/'plots').glob('*_dominant_zoom.png'));assert len(zooms)==len(ids) and {p.name.removesuffix('_dominant_zoom.png') for p in zooms}==set(ids)
 assert all(p.stat().st_size>20000 for p in zooms)
 overview=OUT/'plots/all72_censat_outside_overview.png';assert overview.exists() and overview.stat().st_size>100000
 index=['# 72 k-mer landscape figures','', 'CenSat outside exact-hit数が多い順。個別図にはCHM13座標、annotation、hap group、group × chromosome countを収録。拡大図はCHM13で最もヒットが多い染色体。','', '| k-mer | family | 用途 | CHM13 exact | CenSat外 | 個別図 | 拡大図 |','|---|---|---|---:|---:|---|---|']
 for r in candidates:
  kid=r['kmer_id'];index.append(f"| `{kid}` | {r['target_family']} | {r['candidate_type']} | {int(r['CHM13_total_exact_hits']):,} | {int(r['CHM13_outside_broad_CenSat']):,} | [PNG](plots/{kid}_landscape.png) | [zoom](plots/{kid}_dominant_zoom.png) |")
 (OUT/'PLOT_INDEX.md').write_text('\n'.join(index)+'\n')
 validation={'status':'PASS','STEP4_complete_sha256':source,'source_STEP4_unchanged':sha(ROOT/'step4/COMPLETE.json')==source,'queries':len(ids),'CHM13_exact_hits':s['CHM13_exact_hits'],'hap_chromosome_rows':len(ids)*nh*24,'all_72_per_query_PNG_verified':True,'all_72_dominant_zoom_PNG_verified':True,'CHM13_outside_overview_PNG_verified':True,'original_top3_coordinate_sets_reproduced_by_one_pass_scanner':True,'all_72_fine_HSat2_HSat3_counts_match_STEP4':True,'outside_exact_hits_total':sum(int(x['outside_broad_CenSat']) for x in totals),'queries_with_outside_hits':sum(int(x['outside_broad_CenSat'])>0 for x in totals),'finished_at':time.strftime('%Y-%m-%dT%H:%M:%S%z'),'script_sha256':{p.name:sha(p) for p in [ROOT/'scripts/scan72_chm13.cpp',ROOT/'scripts/all72_landscape.py',ROOT/'scripts/plot_all72_landscape.R',ROOT/'scripts/plot_all72_zoom.R',ROOT/'scripts/finalize_all72_landscape.py']}}
 (OUT/'COMPLETE.json').write_text(json.dumps(validation,indent=2)+'\n');print(json.dumps(validation,indent=2))
if __name__=='__main__':main()
