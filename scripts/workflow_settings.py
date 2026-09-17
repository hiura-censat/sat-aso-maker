"""Optional settings shared by the standalone scripts and the Snakemake run."""
import os
from pathlib import Path


def settings():
    path = os.environ.get('SAT_ASO_CONFIG')
    if not path:
        return {
            'selection': {
                'selected_E': 10,
                'pooled_min_count': 10,
                'hap_min_count': 500,
                'broad_hap_fraction': 0.90,
            },
            'samtools': '/home/senescence/miniconda3/envs/kmer_validation_env/bin/samtools',
            'invariants': {'k': 16, 'max_mismatches': 2},
        }
    import yaml
    data = yaml.safe_load(Path(path).read_text())
    select = data['selection']
    if not isinstance(select['selected_E'], (int, float)) or not 0 < select['selected_E'] < 100:
        raise ValueError('selection.selected_E must be a positive finite threshold below 100')
    if not all(isinstance(select[x], int) for x in ('pooled_min_count', 'hap_min_count')):
        raise ValueError('count minima must be integers')
    if select['pooled_min_count'] < 1 or select['hap_min_count'] < 1:
        raise ValueError('count minima must be positive')
    if not 0 < select['broad_hap_fraction'] <= 1:
        raise ValueError('selection.broad_hap_fraction must be in (0, 1]')
    if data['invariants'] != {'k': 16, 'max_mismatches': 2}:
        raise ValueError('the current scanners require k=16 and max_mismatches=2')
    return data
