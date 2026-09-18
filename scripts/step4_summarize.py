#!/usr/bin/env python3
"""Merge reused STEP3 and new STEP4 mismatch counts and rank 72 ASO targets."""
import csv, gzip, json, math, re
from collections import defaultdict
from pathlib import Path
import numpy as np
from step3_prepare import ROOT
from step2_pipeline import MANIFEST, CHROMS
from workflow_settings import settings

S3=ROOT/'step3'; S4=ROOT/'step4'; M=S4/'mismatch'
REGIONS=['HSat2','HSat3','HSat23','background']

def read(path):
    with (gzip.open(path,'rt') if str(path).endswith('.gz') else open(path)) as f:
        return list(csv.DictReader(f,delimiter='\t'))
def write(path,rows,fields=None):
    rows=list(rows)
    if fields is None: fields=list(rows[0]) if rows else []
    with (gzip.open(path,'wt') if str(path).endswith('.gz') else open(path,'w')) as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',extrasaction='ignore');w.writeheader();w.writerows(rows)
def fnum(x):
    try:return float(x)
    except:return float('nan')
def escore(a,ad,b,bd): return float(np.log2((a/ad+1e-9)/(b/bd+1e-9)))
def rev(s): return s.translate(str.maketrans('ACGT','TGCA'))[::-1]
def hdist(a,b): return sum(x!=y for x,y in zip(a,b))

def load_merged():
    q=read(M/'queries.tsv'); old=read(S3/'candidates/mismatch_queries.tsv'); new=read(M/'new_queries.tsv')
    oldix={r['kmer_id']:i for i,r in enumerate(old)}; newix={r['kmer_id']:i for i,r in enumerate(new)}
    nq=len(q); nh=len(MANIFEST); nc=len(CHROMS)
    counts=np.zeros((nq,nh,nc,4,3),dtype=np.uint64); den=np.zeros((nh,nc,4),dtype=np.uint64)
    pos=np.zeros((nq,nh,4,3,16),dtype=np.uint64); avail=np.zeros((nh,nc),bool)
    for hi,r in enumerate(MANIFEST):
        h=r['hap']; cache={}
        for origin,path in [('STEP3',S3/'mismatch'),('STEP4',M)]:
            qc=json.loads((path/'counts'/(h+'.json')).read_text())
            assert qc['zero_mismatch_matches_step1'] and qc['all_valid_windows_match_step0'] and qc['mismatch_position_accounting']
            z=np.load(path/'counts'/(h+'.npz')); cache[origin]=z
        ch=[str(x) for x in cache['STEP4']['chromosomes']]; ci=[CHROMS.index(x) for x in ch]
        np.testing.assert_array_equal(cache['STEP3']['valid_starts'],cache['STEP4']['valid_starts'])
        np.testing.assert_array_equal(cache['STEP3']['chromosomes'],cache['STEP4']['chromosomes'])
        den[hi,ci]=cache['STEP4']['valid_starts']; avail[hi,ci]=True
        for qi,x in enumerate(q):
            origin=x['query_origin']; si=(oldix if origin=='STEP3' else newix)[x['kmer_id']]; z=cache[origin]
            counts[qi,hi,ci]=z['counts'][si]; pos[qi,hi]=z['position_counts'][si]
        for z in cache.values(): z.close()
    return q,counts,den,pos,avail

def main():
    q,c,den,pos,avail=load_merged(); nq,nh,nc,_,_=c.shape
    ref=np.array([r['sample']=='CHM13' for r in MANIFEST]); non=~ref
    broad_min=int(np.ceil(settings()['selection']['broad_hap_fraction']*int(non.sum())))
    total=c.sum(axis=2); pool=c[:,non].sum(axis=(1,2)); poolden=den[non].sum(axis=(0,1)); cum=pool.cumsum(axis=2); local=total.cumsum(axis=3)
    pooled=[]; robust=[]
    for qi,x in enumerate(q):
        fam=x['target_family']; fi=['HSat2','HSat3'].index(fam); oi=1-fi; flags=[]
        for radius in range(3):
            t,o,b=map(int,[cum[qi,fi,radius],cum[qi,oi,radius],cum[qi,3,radius]])
            eb=escore(t,poolden[fi],b,poolden[3]); ef=escore(t,poolden[fi],o,poolden[oi]); frac=t/(t+o) if t+o else float('nan')
            with np.errstate(divide='ignore',invalid='ignore'):
                le=np.log2((local[qi,:,fi,radius]/den.sum(axis=1)[:,fi]+1e-9)/(local[qi,:,3,radius]/den.sum(axis=1)[:,3]+1e-9))
                lf=np.log2((local[qi,:,fi,radius]/den.sum(axis=1)[:,fi]+1e-9)/(local[qi,:,oi,radius]/den.sum(axis=1)[:,oi]+1e-9))
                lfrac=local[qi,:,fi,radius]/(local[qi,:,fi,radius]+local[qi,:,oi,radius])
            good=non&(local[qi,:,fi,radius]>=10)&(le>10)&(lf>10)&(lfrac>=.999)
            broad=int(good.sum()); passed=eb>10 and ef>10 and frac>=.999 and broad>=broad_min; flags.append(passed)
            pooled.append(dict(query_index=qi,kmer_id=x['kmer_id'],canonical_kmer=x['canonical_kmer'],target_family=fam,query_category=x['query_category'],max_mismatches=radius,target_count_le_radius=t,other_family_count_le_radius=o,background_count_le_radius=b,target_vs_background_E=eb,target_vs_other_family_E=ef,target_family_count_fraction=frac,strict_passing_nonreference_haps=broad,passes_strict_pan_criteria=int(passed)))
        tier='A_le2mm' if flags[2] else ('B_le1mm' if flags[1] else ('C_exact' if flags[0] else 'D_fails_exact_strict'))
        robust.append(dict(kmer_id=x['kmer_id'],strict_pan_pass_0mm=int(flags[0]),strict_pan_pass_le1mm=int(flags[1]),strict_pan_pass_le2mm=int(flags[2]),robustness_tier=tier))
    write(M/'pooled_mismatch_scores.tsv',pooled); write(M/'candidate_robustness.tsv',robust)

    # Conditional chromosome concentration within available annotated target-family regions.
    chromosome_rows=[]; cc=c[:,non].sum(axis=1).cumsum(axis=3); lc=c.cumsum(axis=4)
    for qi,x in enumerate(q):
        fi=['HSat2','HSat3'].index(x['target_family'])
        for radius in range(3):
            v=cc[qi,:,fi,radius]; tot=int(v.sum()); best=int(v.argmax()); tied=(v==v[best]).sum()>1
            per=lc[qi,:,:,fi,radius]; mx=per.max(axis=1); bi=per.argmax(axis=1); lt=(per==mx[:,None]).sum(axis=1)>1
            measured=den[:,:,fi]>0; eligible=non&(per.sum(axis=1)>0)&measured[:,best]&(measured.sum(axis=1)>=2)
            if tied:eligible[:]=False
            support=int((eligible&~lt&(bi==best)).sum()); n=int(eligible.sum()); frac=float(v[best]/tot) if tot else float('nan'); rec=support/n if n else float('nan')
            chromosome_rows.append(dict(query_index=qi,kmer_id=x['kmer_id'],target_family=x['target_family'],query_category=x['query_category'],max_mismatches=radius,dominant_chromosome=CHROMS[best] if tot else 'NA',max_target_count_fraction=frac,informative_haps=n,same_unique_dominant_chromosome_haps=support,recurrence_fraction=rec,passes_chromosome_screen=int(frac>=.9 and rec>=.8 and n>=100)))
    write(M/'chromosome_robustness.tsv',chromosome_rows)

    # Group contrasts use only the chromosome on which the group was fitted.
    groups=[]; summaries=json.loads((S3/'clustering_summary.json').read_text())
    for case in summaries:
        if not case['supported_partition']:continue
        fam=case['family'];case_chromosome=case['chromosome'];fi=['HSat2','HSat3'].index(fam);ci=CHROMS.index(case_chromosome)
        with np.load(S3/'clustering'/case['case']/'model.npz') as z: hi=z['hap_indices'];labels=z['labels']
        assert np.all(den[hi,ci,fi]>0) and not ref[hi].any()
        for qi,x in enumerate(q):
            if x['target_family']!=fam:continue
            ct=c[qi,hi,ci,fi,:].cumsum(axis=1);density=ct/den[hi,ci,fi,None]*1e6
            for g in np.unique(labels):
                inside=labels==g
                for radius in range(3):
                    pin=float((ct[inside,radius]>=10).mean());pout=float((ct[~inside,radius]>=10).mean())
                    aa=float(np.median(density[inside,radius]));bb=float(np.median(density[~inside,radius]))
                    ratio=float(np.log2((aa+.001)/(bb+.001)))
                    groups.append(dict(query_index=qi,kmer_id=x['kmer_id'],target_family=fam,query_category=x['query_category'],chromosome=case_chromosome,group=f'G{g}',scope='chromosome',max_mismatches=radius,in_group_count_ge10_fraction=pin,out_group_count_ge10_fraction=pout,in_group_median_density_per_M=aa,out_group_median_density_per_M=bb,log2_median_density_ratio=ratio,passes_group_enrichment_screen=int(pin>=.8 and ratio>=2),passes_group_specific_screen=int(pin>=.8 and pout<=.2 and ratio>=2)))
    write(M/'group_robustness.tsv',groups)

    # Similarity components: Hamming <=2 in either orientation or any shared 12-mer.
    parent=list(range(nq))
    def find(a):
        while parent[a]!=a: parent[a]=parent[parent[a]];a=parent[a]
        return a
    def union(a,b):
        a,b=find(a),find(b)
        if a!=b:parent[b]=a
    seqs=[x['canonical_kmer'] for x in q]
    for i in range(nq):
        for j in range(i):
            s,t=seqs[i],seqs[j]; rt=rev(t)
            if min(hdist(s,t),hdist(s,rt))<=2 or any(s[k:k+12] in t or s[k:k+12] in rt for k in range(5)): union(i,j)
    roots={}; classes=[]
    for i in range(nq): roots.setdefault(find(i),len(roots)+1); classes.append(f'SC{roots[find(i)]:03d}')

    exact={r['kmer_id']:r for r in read(S4/'ranked/all_exact_ranked.tsv.gz')}; pidx={(r['kmer_id'],int(r['max_mismatches'])):r for r in pooled}; ridx={r['kmer_id']:r for r in robust}; chidx={(r['kmer_id'],int(r['max_mismatches'])):r for r in chromosome_rows}; gidx={(r['kmer_id'],r['chromosome'],r['group'],int(r['max_mismatches'])):r for r in groups}
    ranked=[]
    for qi,x in enumerate(q):
        e=exact[x['kmer_id']]; p2=pidx[(x['kmer_id'],2)]; typ='pan' if x['query_category'].startswith('pan') else ('chromosome' if x['query_category'].startswith('chromosome') else 'regional_group')
        category_pass=[bool(ridx[x['kmer_id']][f'strict_pan_pass_{"0mm" if d==0 else "le"+str(d)+"mm"}']) for d in range(3)]; category_note='strict pan criteria'
        if typ=='chromosome':
            m=re.search(r'_chr(?:omosome_)?(\w+)$',x['query_category']); wanted='chr'+m.group(1) if m and not m.group(1).startswith('chr') else (m.group(1) if m else '')
            category_pass=[bool(chidx[(x['kmer_id'],d)]['passes_chromosome_screen']) and (not wanted or chidx[(x['kmer_id'],d)]['dominant_chromosome']==wanted) for d in range(3)];category_note=f'expected {wanted}' if wanted else 'dominant chromosome from mismatch count'
        elif typ=='regional_group':
            m=re.search(r'(chr(?:[0-9]+|X|Y))_(G[0-9]+)',x['query_category'])
            if not m:raise ValueError('regional category lacks chromosome and group: '+x['query_category'])
            wanted_chromosome,wanted=m.groups();specific='specific' in x['query_category']
            category_pass=[bool(gidx[(x['kmer_id'],wanted_chromosome,wanted,d)]['passes_group_specific_screen' if specific else 'passes_group_enrichment_screen']) for d in range(3)]
            category_note=f'{wanted_chromosome} {wanted} '+('specific' if specific else 'enriched')
        cat_tier='A_le2mm' if category_pass[2] else ('B_le1mm' if category_pass[1] else ('C_exact' if category_pass[0] else 'D_fails_exact'))
        # Smooth mismatch score, bounded and transparent; strict tiers remain separate columns.
        mm=(min(1,max(0,fnum(p2['target_vs_background_E'])/10))*.30 + min(1,max(0,fnum(p2['target_vs_other_family_E'])/10))*.25 + min(1,max(0,(fnum(p2['target_family_count_fraction'])-.9)/.1))*.20 + min(1,int(p2['strict_passing_nonreference_haps'])/broad_min)*.25)
        score=.55*fnum(e['exact_composite_score'])+.35*mm+.10*(1 if category_pass[2] else (.6 if category_pass[1] else (.3 if category_pass[0] else 0)))
        ranked.append(dict(kmer_id=x['kmer_id'],canonical_target_5to3=x['canonical_kmer'],reverse_complement_target_5to3=x['reverse_complement'],aso_5to3_if_canonical_target_is_transcribed=x['reverse_complement'],aso_5to3_if_reverse_complement_target_is_transcribed=x['canonical_kmer'],target_family=x['target_family'],candidate_type=typ,query_category=x['query_category'],query_origin=x['query_origin'],sequence_class=classes[qi],final_score=score,exact_composite_score=e['exact_composite_score'],mismatch_score=mm,strict_pan_robustness_tier=ridx[x['kmer_id']]['robustness_tier'],category_robustness_tier=cat_tier,category_interpretation=category_note,passes_category_0mm=int(category_pass[0]),passes_category_le1mm=int(category_pass[1]),passes_category_le2mm=int(category_pass[2]),target_vs_background_E_le2mm=p2['target_vs_background_E'],target_vs_other_family_E_le2mm=p2['target_vs_other_family_E'],target_family_fraction_le2mm=p2['target_family_count_fraction'],strict_passing_haps_le2mm=p2['strict_passing_nonreference_haps'],GC_percent=e['GC_percent'],max_homopolymer=e['max_homopolymer'],self_complementary_stem_proxy=e['self_complementary_stem_proxy']))
    ranked.sort(key=lambda r:(r['target_family'],r['candidate_type'],-r['final_score']))
    seen=defaultdict(int)
    for r in ranked: seen[(r['target_family'],r['candidate_type'])]+=1;r['rank_within_family_and_type']=seen[(r['target_family'],r['candidate_type'])]
    fields=list(ranked[0]);write(S4/'ranked/mismatch_evaluated_ranked.tsv',ranked,fields)
    for fam in ['HSat2','HSat3']:
        for typ in ['pan','regional_group','chromosome']:
            rr=[r for r in ranked if r['target_family']==fam and r['candidate_type']==typ]
            write(S4/'shortlists'/f'{fam}_{typ}_candidates.tsv',rr[:20],fields)
    pp=pos[:,non].sum(axis=1)
    write(M/'mismatch_position_counts.tsv',(dict(query_index=qi,kmer_id=x['kmer_id'],region=REGIONS[ri],exact_distance=d,query_position_1based=k+1,weighted_count=int(pp[qi,ri,d,k])) for qi,x in enumerate(q) for ri in range(4) for d in [1,2] for k in range(16)))
    for d in range(3): np.testing.assert_array_equal(pp[:,:,d,:].sum(axis=2),pool[:,:,d]*d)
    np.savez_compressed(M/'combined_hap_counts.npz',counts=total,valid_starts=den.sum(axis=1),chromosome_counts=c,chromosome_valid_starts=den,chromosome_available=avail,position_counts=pos,query_ids=np.array([x['kmer_id'] for x in q]),haplotypes=np.array([r['hap'] for r in MANIFEST]),is_reference=ref,chromosomes=np.array(CHROMS))
    summary={'status':'PASS','queries':nq,'reused_STEP3_queries':sum(x['query_origin']=='STEP3' for x in q),'new_STEP4_queries':sum(x['query_origin']=='STEP4' for x in q),'sequence_sets':nh,'nonreference_haps':int(non.sum()),'broad_haps_min':broad_min,'strict_pan_pass_counts_by_radius':[sum(int(r[f'strict_pan_pass_{"0mm" if d==0 else "le"+str(d)+"mm"}']) for r in robust) for d in range(3)],'category_pass_counts_by_radius':[sum(int(r[f'passes_category_{"0mm" if d==0 else "le"+str(d)+"mm"}']) for r in ranked) for d in range(3)],'robustness_tiers':dict((t,sum(r['category_robustness_tier']==t for r in ranked)) for t in ['A_le2mm','B_le1mm','C_exact','D_fails_exact']),'method':'exhaustive Hamming radius <=2 against both orientations; no indels','scope':'eligible T2T chromosomes and existing HSat2/3 BED annotations'}
    (M/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))
if __name__=='__main__': main()
