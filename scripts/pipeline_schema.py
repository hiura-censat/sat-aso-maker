"""Shared tabular schemas and fail-fast validation for workflow stages."""

STEP2_CANDIDATE_COLUMNS = [
    'row_index', 'kmer_id', 'canonical_kmer', 'reverse_complement',
    'pooled_selected_threshold', 'pooled_and_broad_threshold',
]

STEP2_METRIC_COLUMNS = [
    'kmer_id', 'canonical_kmer', 'pooled_selected_threshold',
    'pooled_and_broad_threshold', 'HSat2_count', 'HSat3_count',
    'HSat23_union_count', 'HSat2_count_fraction', 'HSat2_density_per_M',
    'HSat3_density_per_M', 'log2_HSat2_vs_HSat3_density',
    'HSat2_prevalence_haps', 'HSat3_prevalence_haps',
    'union_prevalence_haps', 'union_prevalence_fraction',
    'union_count_ge10_haps', 'union_count_ge100_haps', 'target_count_p10',
    'target_count_median', 'target_count_p90', 'target_density_per_M_p10',
    'target_density_per_M_median', 'target_density_per_M_p90',
    'target_density_CV', 'family_informative_haps',
    'HSat2_density_favored_haps', 'HSat3_density_favored_haps',
]

def require_columns(actual, required, source):
    actual = list(actual or [])
    missing = [name for name in required if name not in actual]
    if missing:
        raise ValueError(f'{source} schema mismatch; missing columns: {missing}; actual: {actual}')
    return actual
