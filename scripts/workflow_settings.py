"""Optional settings shared by the standalone scripts and the Snakemake run."""
import os
import math
from pathlib import Path


def settings():
    path = os.environ.get('SAT_ASO_CONFIG')
    if not path:
        return {
            'selection': {
                'selected_E': 10,
                'sweep_E': [2, 4, 6, 8, 10, 12, 14, 16, 18, 20],
                'pooled_min_count': 100,
                'hap_min_count': 10,
                'broad_hap_fraction': 0.90,
            },
            'samtools': '/home/senescence/miniconda3/envs/kmer_validation_env/bin/samtools',
            'invariants': {'k': 16, 'max_mismatches': 2},
        }
    import yaml
    data = yaml.safe_load(Path(path).read_text())
    select = data['selection']
    thresholds = select['sweep_E']
    if not thresholds or not all(isinstance(x, (int, float)) and math.isfinite(x) for x in thresholds) or thresholds != sorted(set(thresholds)):
        raise ValueError('selection.sweep_E must be a sorted unique nonempty list')
    if select['selected_E'] not in thresholds:
        raise ValueError('selection.selected_E must be in selection.sweep_E')
    if not all(isinstance(select[x], int) for x in ('pooled_min_count', 'hap_min_count')):
        raise ValueError('count minima must be integers')
    if select['pooled_min_count'] < 100 or select['hap_min_count'] < 10:
        raise ValueError('count minima below 100/10 require changing the STEP1 preliminary universe')
    if not 0 < select['broad_hap_fraction'] <= 1:
        raise ValueError('selection.broad_hap_fraction must be in (0, 1]')
    if data['invariants'] != {'k': 16, 'max_mismatches': 2}:
        raise ValueError('the current scanners require k=16 and max_mismatches=2')
    return data
