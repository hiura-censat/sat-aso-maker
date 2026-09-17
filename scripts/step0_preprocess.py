#!/usr/bin/env python3
"""Build HSat2/3 region views and count valid 16-mer windows without copying genomes."""
import argparse, csv, gzip, json, re, shutil, os, subprocess
from pathlib import Path
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed

ROOT=Path(__file__).resolve().parents[1]
SOURCE=ROOT.parent/'VallePrep_v0.0.0/results/data'
OUT=ROOT/'step0'

def write(path, header, rows):
    path.parent.mkdir(parents=True,exist_ok=True)
    with path.open('w') as f:
        w=csv.writer(f,delimiter='\t'); w.writerow(header); w.writerows(rows)

def merge(xs):
    result=[]
    for s,e in sorted(xs):
        if result and s<=result[-1][1]: result[-1]=(result[-1][0],max(e,result[-1][1]))
        else: result.append((s,e))
    return result

def intersect(a,b):
    out=[]; i=j=0
    while i<len(a) and j<len(b):
        s=max(a[i][0],b[j][0]); e=min(a[i][1],b[j][1])
        if s<e: out.append((s,e))
        if a[i][1]<b[j][1]: i+=1
        else: j+=1
    return out

def subtract(a,b):
    out=[]
    for s,e in a:
        for x,y in b:
            if y<=s: continue
            if x>=e: break
            if s<x: out.append((s,x))
            s=max(s,y)
            if s>=e: break
        if s<e: out.append((s,e))
    return out

def records(path):
    with gzip.open(str(path),'rb') as f:
        name=None; chunks=[]
        for line in f:
            if line.startswith(b'>'):
                if name is not None: yield name,b''.join(chunks)
                name=line[1:].split()[0].decode(); chunks=[]
            else: chunks.append(line.strip())
        if name is not None: yield name,b''.join(chunks)

def stats(job):
    hap,fa,lengths,regions=job
    dest=OUT/'qc'/f'{hap}.stats.tsv'
    if dest.exists(): return hap,'cached'
    scanner=ROOT/'scripts/step0_scan'
    region_dir=OUT/'regions'/hap
    if scanner.exists() and all((region_dir/(f+'.bed')).exists() for f in ['HSat2','HSat3','HSat23','ambiguous','background']):
        dest.parent.mkdir(parents=True,exist_ok=True)
        tmp=dest.with_suffix('.tmp')
        subprocess.run([str(scanner),fa,str(region_dir),hap,str(tmp)],check=True)
        tmp.replace(dest)
        return hap,'done'
    rows=[]; seen=set()
    for chrom,seq in records(fa):
        if chrom not in lengths: continue
        if chrom in seen: raise ValueError('duplicate FASTA chromosome '+hap+':'+chrom)
        seen.add(chrom)
        if len(seq)!=lengths[chrom]: raise ValueError('FASTA/FAI length mismatch '+hap+':'+chrom)
        for family,bychr in regions.items():
            intervals=bychr.get(chrom,[]); valid=0; acgt=0
            for s,e in intervals:
                for match in re.finditer(b'[ACGTacgt]+',seq[s:e]):
                    n=match.end()-match.start(); acgt+=n; valid+=max(0,n-15)
            total=sum(e-s for s,e in intervals)
            rows.append([hap,chrom,family,len(intervals),total,acgt,total-acgt,valid])
    if seen!=set(lengths): raise ValueError('missing FASTA chromosomes '+hap)
    tmp=dest.with_suffix('.tmp'); write(tmp,STAT_HEADER,rows); tmp.replace(dest)
    return hap,'done'

STAT_HEADER=['hap','chromosome','region','interval_count','bp','ACGT_bp','non_ACGT_bp','valid_16mer_starts']

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--workers',type=int,default=4); args=ap.parse_args()
    jobs=[]; manifest=[]; statuses=[]; qc=[]; invalid=[]
    for passbed in sorted(SOURCE.glob('*/centro/*.t2t.pass.bed')):
        hap=passbed.name[:-len('.t2t.pass.bed')]; sample=passbed.parent.parent.name; base=passbed.parent.parent
        fa=base/'filfa'/f'{hap}.t2t.filtered.fasta.gz'; fai=Path(str(fa)+'.fai')
        bed=base/'hsat23'/f'{hap}.HSat2and3_Regions.lifted.bed'; log=base/'hsat23'/f'{hap}.HSat2and3_Regions.HSat2and3.log'
        passed={l.split()[0] for l in passbed.read_text().splitlines() if l.strip() and not l.startswith('#')}
        lengths={l.split()[0]:int(l.split()[1]) for l in fai.read_text().splitlines()} if fai.exists() else {}
        text=log.read_text(errors='replace') if log.exists() else ''
        reason='PASS' if all(p.exists() for p in [fa,fai,bed,log]) and 'printed output to ' in text and 'Total running time:' in text and not re.search(r'\bERROR\b',text) else 'EXCLUDE_annotation_incomplete_or_missing_input'
        if not passed<=set(lengths): reason='EXCLUDE_missing_filtered_chromosome'
        ann=defaultdict(lambda:defaultdict(list)); rawcount=0
        if bed.exists():
            for n,line in enumerate(bed.read_text().splitlines(),1):
                if not line.strip() or line.startswith('#'): continue
                f=line.split(); rawcount+=1
                if len(f)<6 or f[3] not in ('HSat2','HSat3') or f[5] not in ('+','-'): raise ValueError(f'{bed}:{n}: bad annotation')
                c,s,e=f[0],int(f[1]),int(f[2])
                if c not in lengths or not 0<=s<e<=lengths[c]:
                    invalid.append([hap,str(bed),n,c,s,e]); reason='EXCLUDE_invalid_annotation_coordinate'; continue
                ann[f[3]][c].append((s,e))
        manifest.append([sample,hap,reason,str(fa),str(fai),str(passbed),str(bed),str(log),len(passed),rawcount])
        for c in sorted(set(lengths)|passed):
            statuses.append([sample,hap,c,int(c in passed),int(c in lengths),reason if c in passed else 'EXCLUDE_not_t2t_pass',len(ann['HSat2'][c]),len(ann['HSat3'][c])])
        if reason!='PASS': continue
        lengths={c:lengths[c] for c in sorted(passed)}
        regions={k:{} for k in ['HSat2','HSat3','HSat23','ambiguous','background']}
        for c,L in lengths.items():
            a=merge(ann['HSat2'][c]); b=merge(ann['HSat3'][c]); overlap=intersect(a,b); union=merge(a+b)
            regions['HSat2'][c]=subtract(a,overlap); regions['HSat3'][c]=subtract(b,overlap)
            regions['HSat23'][c]=union; regions['ambiguous'][c]=overlap; regions['background'][c]=subtract([(0,L)],union)
            assert sum(e-s for s,e in union+regions['background'][c])==L
        folder=OUT/'regions'/hap; folder.mkdir(parents=True,exist_ok=True)
        shutil.copyfile(bed,folder/'source_annotation.bed')
        for family,bychr in regions.items():
            with (folder/f'{family}.bed').open('w') as f:
                for c,xs in bychr.items():
                    for s,e in xs: f.write(f'{c}\t{s}\t{e}\t{family}\n')
        seqdir=OUT/'sequences'/hap; seqdir.mkdir(parents=True,exist_ok=True)
        links=[(fa,seqdir/'source.fasta.gz'),(fai,seqdir/'source.fasta.gz.fai')]
        gzi=Path(str(fa)+'.gzi')
        if gzi.exists(): links.append((gzi,seqdir/'source.fasta.gz.gzi'))
        for src,dst in links:
            if not dst.is_symlink(): dst.symlink_to(src)
        qc.append([hap,'PASS',rawcount,sum(len(v) for v in regions['ambiguous'].values()),sum(e-s for xs in regions['ambiguous'].values() for s,e in xs)])
        jobs.append((hap,str(fa),lengths,regions))
    write(OUT/'invalid_annotations.tsv',['hap','bed','line','chromosome','start','end'],invalid)
    write(OUT/'manifest.tsv',['sample','hap','status','fasta','fai','t2t_pass_bed','hsat23_lifted_bed','altemose_log','t2t_chromosomes','annotation_records'],manifest)
    write(OUT/'chromosome_status.tsv',['sample','hap','chromosome','t2t_pass','in_filtered_fasta','status','HSat2_raw_regions','HSat3_raw_regions'],statuses)
    write(OUT/'qc_summary.tsv',['hap','coordinate_status','source_records','ambiguous_intervals','ambiguous_bp'],qc)
    metadata={}
    for path in sorted({value for row in manifest for value in row[3:8]}):
        p=Path(path)
        metadata[path]={'size':p.stat().st_size,'mtime_ns':p.stat().st_mtime_ns} if p.exists() else None
    (OUT/'input_metadata.json').write_text(json.dumps(metadata,indent=2)+'\n')
    views=[]
    for hap,fa,_,regions in jobs:
        for family in regions:
            views.append([hap,family,fa,str(OUT/'regions'/hap/(family+'.bed')),'filtered_oriented_fasta_0based_halfopen'])
    write(OUT/'sequence_views.tsv',['hap','region','source_fasta','region_bed','coordinate_system'],views)
    print(f'Eligible {len(jobs)} / {len(manifest)}; scanning FASTA with {args.workers} workers',flush=True)
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        for i,future in enumerate(as_completed([pool.submit(stats,j) for j in jobs]),1):
            hap,status=future.result(); print(f'{i}/{len(jobs)} {hap} {status}',flush=True)
    with (OUT/'region_stats.tsv').open('w') as f:
        f.write('\t'.join(STAT_HEADER)+'\n')
        for hap,*_ in jobs:
            with (OUT/'qc'/f'{hap}.stats.tsv').open() as source:
                next(source); shutil.copyfileobj(source,f)
    (OUT/'COMPLETE.json').write_text(json.dumps({'haps':len(jobs),'excluded':len(manifest)-len(jobs),'k':16,'sequence_storage':'source FASTA symlinks plus region BED views'},indent=2)+'\n')

if __name__=='__main__': main()
