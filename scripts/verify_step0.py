#!/usr/bin/env python3
import csv,json
from collections import Counter,defaultdict
from pathlib import Path
from step0_preprocess import OUT

def read(path):
    with path.open() as f:return list(csv.DictReader(f,delimiter='\t'))
manifest=read(OUT/'manifest.tsv'); eligible={r['hap']:r for r in manifest if r['status']=='PASS'}
rows=read(OUT/'region_stats.tsv'); groups=defaultdict(dict); totals=Counter()
for r in rows:
    h,c,f=r['hap'],r['chromosome'],r['region']
    assert h in eligible and f not in groups[h,c],(h,c,f)
    n={k:int(r[k]) for k in ['interval_count','bp','ACGT_bp','non_ACGT_bp','valid_16mer_starts']}
    assert all(v>=0 for v in n.values()) and n['bp']==n['ACGT_bp']+n['non_ACGT_bp']
    assert n['valid_16mer_starts']<=n['ACGT_bp']
    groups[h,c][f]=n
    totals[f+'_bp']+=n['bp'];totals[f+'_valid_16mer_starts']+=n['valid_16mer_starts']
expected={(r['hap'],r['chromosome']) for r in read(OUT/'chromosome_status.tsv') if r['status']=='PASS'}
assert set(groups)==expected
lengths={}
for h,r in eligible.items():
    with open(r['fai']) as f:lengths[h]={a[0]:int(a[1]) for a in (l.split() for l in f)}
for (h,c),g in groups.items():
    assert set(g)=={'HSat2','HSat3','HSat23','ambiguous','background'}
    assert g['HSat23']['bp']+g['background']['bp']==lengths[h][c]
    for field in ['bp','ACGT_bp','non_ACGT_bp']:
        assert g['HSat23'][field]==sum(g[f][field] for f in ['HSat2','HSat3','ambiguous'])
for path,old in json.loads((OUT/'input_metadata.json').read_text()).items():
    p=Path(path)
    if old is None: assert not p.exists(),path
    else: assert {'size':p.stat().st_size,'mtime_ns':p.stat().st_mtime_ns}==old,path
summary={'eligible_sequence_sets':len(eligible),'includes_CHM13_reference':any(r['sample']=='CHM13' for r in eligible.values()),'samples':len({r['sample'] for r in eligible.values()}),'hap_chromosome_pairs':len(groups),'statistics_rows':len(rows),'excluded':[{'hap':r['hap'],'reason':r['status']} for r in manifest if r['status']!='PASS'],'totals':dict(totals),'input_metadata_unchanged':True}
(OUT/'validation_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps(summary,indent=2))
