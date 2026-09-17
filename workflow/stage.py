#!/usr/bin/env python3
"""Run one existing pipeline stage inside a Snakemake work directory."""
import argparse
import csv
import json
import os
import subprocess
from pathlib import Path

import yaml


def manifest_haps(run):
    with (run / 'step0/manifest.tsv').open() as stream:
        return [row['hap'] for row in csv.DictReader(stream, delimiter='\t') if row['status'] == 'PASS']


def check_shards(run, folder, suffix):
    missing = [h for h in manifest_haps(run) if not (run / folder / (h + suffix)).is_file()]
    if missing:
        raise RuntimeError(f'{len(missing)} missing {folder} shards, e.g. {missing[:3]}')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('stage')
    parser.add_argument('--run', required=True)
    parser.add_argument('--threads', type=int, default=1)
    args = parser.parse_args()
    run = Path(args.run).resolve()
    cfg = yaml.safe_load((run / 'workflow_config.yaml').read_text())
    py = cfg['python']
    env = os.environ.copy()
    env['SAT_ASO_CONFIG'] = str(run / 'workflow_config.yaml')
    env['SAT_ASO_PYTHON'] = py
    env['SAT_ASO_STEP3_WORKERS'] = str(args.threads)
    env['SAT_ASO_SELECTED_E'] = str(cfg['selection']['selected_E'])
    env['SAT_ASO_POOLED_MIN'] = str(cfg['selection']['pooled_min_count'])
    env['SAT_ASO_HAP_MIN'] = str(cfg['selection']['hap_min_count'])
    env['SAT_ASO_BROAD_FRACTION'] = str(cfg['selection']['broad_hap_fraction'])
    env['OPENBLAS_NUM_THREADS'] = '2'
    env['OMP_NUM_THREADS'] = '2'
    script = lambda name: [py, 'scripts/' + name]
    commands = {
        'step0': [script('test_step0.py'), script('step0_preprocess.py') + ['--workers', str(args.threads)], script('verify_step0.py')],
        'step1_discovery': [script('test_step1.py'), script('step1_pipeline.py') + ['--phase', 'discovery', '--workers', str(args.threads)]],
        'step1_counts': [script('step1_pipeline.py') + ['--phase', 'counts', '--workers', str(args.threads)]],
        'step1_finish': [script('finish_step1.py')],
        'sweep': [script('sweep_step1.py'), ['Rscript', 'scripts/plot_step1_sweep.R']],
        'step2_counts': [script('test_step2.py'), script('step2_pipeline.py') + ['--workers', str(args.threads)]],
        'step2_finish': [script('finish_step2.py')],
        'step3_prepare': [script('step3_prepare.py'), script('test_step3.py')],
        'step3_run': [script('run_step3.py')],
        'step4_rank': [script('step4_rank.py')],
        'step4_neighbors': [script('step4_neighbors.py')],
        'step4_mismatch': [script('step4_mismatch.py') + ['--workers', str(args.threads)]],
        'step4_finish': [script('step4_finalize.py')],
        'top3': [script('top_kmer_landscape.py'), ['Rscript', 'scripts/plot_top_kmer_landscape.R'], script('annotate_top_kmer_hits.py'), ['Rscript', 'scripts/plot_top_kmer_censat_overlay.R']],
        'scan72': [['g++', '-O3', '-std=c++17', 'scripts/scan72_chm13.cpp', '-lz', '-o', 'scripts/scan72_chm13'], ['scripts/scan72_chm13', 'step4/mismatch/queries.tsv', 'step0/sequences/chm13v2.0/source.fasta.gz', 'step4/all72_landscape/raw/chm13_all72_exact_hits.tsv.gz']],
        'landscape': [script('all72_landscape.py')],
        'landscape_plots': [['Rscript', 'scripts/plot_all72_landscape.R'], ['Rscript', 'scripts/plot_all72_zoom.R']],
        'landscape_finish': [script('finalize_all72_landscape.py')],
    }
    if args.stage not in commands:
        raise ValueError(args.stage)
    if args.stage == 'scan72':
        (run / 'step4/all72_landscape/raw').mkdir(parents=True, exist_ok=True)
    for command in commands[args.stage]:
        print('RUN', ' '.join(command), flush=True)
        subprocess.run(command, cwd=run, env=env, check=True)
    if args.stage == 'step1_discovery':
        check_shards(run, 'step1/discovery', '.bin.gz')
    if args.stage == 'step1_counts':
        check_shards(run, 'step1/counts', '.npz')
    if args.stage == 'step2_counts':
        check_shards(run, 'step2/counts', '.npz')
    if args.stage == 'step4_mismatch':
        check_shards(run, 'step4/mismatch/counts', '.npz')
    expected = {
        'step0': ['step0/COMPLETE.json', 'step0/manifest.tsv', 'step0/region_stats.tsv'],
        'step1_discovery': ['step1/preliminary_candidates.u32'],
        'step1_finish': ['step1/COMPLETE.json', 'step1/pooled_counts_and_scores.npz'],
        'sweep': [f"step1/threshold_sweep/candidates_E_gt_{cfg['selection']['selected_E']}.tsv", 'step1/threshold_sweep/validation.json', 'step1/threshold_sweep/threshold_sweep.png'],
        'step2_finish': ['step2/COMPLETE.json', 'step2/matrices/matrix_A_counts.npy'],
        'step3_prepare': ['step3/data/PREPARED.json'],
        'step3_run': ['step3/COMPLETE.json', 'step3/hap_groups.tsv'],
        'step4_rank': ['step4/mismatch/queries.tsv', 'step4/ranked/all_exact_ranked.tsv.gz'],
        'step4_neighbors': ['step4/mismatch/variant_mapping.npz'],
        'step4_finish': ['step4/COMPLETE.json', 'step4/ranked/mismatch_evaluated_ranked.tsv'],
        'top3': ['step4/top_kmer_landscape/censat_overlay_summary.json'],
        'scan72': ['step4/all72_landscape/raw/chm13_all72_exact_hits.tsv.gz'],
        'landscape': ['step4/all72_landscape/summary.json', 'step4/all72_landscape/chm13_candidate_censat_summary.tsv'],
        'landscape_plots': ['step4/all72_landscape/plots/all72_censat_outside_overview.png'],
        'landscape_finish': ['step4/all72_landscape/COMPLETE.json', 'step4/all72_landscape/PLOT_INDEX.md'],
    }
    missing = [name for name in expected.get(args.stage, []) if not (run / name).is_file()]
    if missing:
        raise RuntimeError(f'{args.stage} completed without outputs: {missing}')
    for name in expected.get(args.stage, []):
        if name.endswith('COMPLETE.json') or name.endswith('validation.json'):
            result = json.loads((run / name).read_text())
            if result.get('status') != 'PASS' and name != 'step0/COMPLETE.json':
                raise RuntimeError(f'{name} did not report PASS')
    marker = run / '.workflow' / (args.stage + '.done.json')
    marker.parent.mkdir(exist_ok=True)
    marker.write_text(json.dumps({'stage': args.stage, 'status': 'PASS', 'threads': args.threads}, indent=2) + '\n')


if __name__ == '__main__':
    main()
