import ast
import unittest
from pathlib import Path
from pipeline_schema import STEP2_CANDIDATE_COLUMNS, STEP2_METRIC_COLUMNS, require_columns

ROOT = Path(__file__).resolve().parents[1]

class WorkflowContracts(unittest.TestCase):
    def test_step2_schema_is_complete_and_unique(self):
        self.assertEqual(len(STEP2_CANDIDATE_COLUMNS), len(set(STEP2_CANDIDATE_COLUMNS)))
        self.assertEqual(len(STEP2_METRIC_COLUMNS), len(set(STEP2_METRIC_COLUMNS)))
        require_columns(STEP2_CANDIDATE_COLUMNS, STEP2_CANDIDATE_COLUMNS, 'synthetic candidates')
        require_columns(STEP2_METRIC_COLUMNS, STEP2_METRIC_COLUMNS, 'synthetic metrics')

    def test_active_workflow_has_no_retired_step2_columns(self):
        active = ['step2_pipeline.py', 'summarize_step2.py', 'verify_step2.py',
                  'step3_candidates.py', 'step4_rank.py']
        retired = {'pooled_E_gt_10', 'pooled_and_broad90_E_gt_10'}
        for name in active:
            tree = ast.parse((ROOT / 'scripts' / name).read_text(), filename=name)
            strings = {node.value for node in ast.walk(tree)
                       if isinstance(node, ast.Constant) and isinstance(node.value, str)}
            self.assertFalse(strings & retired, f'{name} uses retired columns {strings & retired}')

    def test_downstream_metric_requirements_exist(self):
        required = {'kmer_id', 'canonical_kmer', 'HSat2_count', 'HSat3_count',
                    'HSat2_count_fraction', 'HSat2_density_per_M',
                    'HSat3_density_per_M', 'log2_HSat2_vs_HSat3_density',
                    'HSat2_prevalence_haps', 'HSat3_prevalence_haps',
                    'union_prevalence_haps'}
        self.assertTrue(required <= set(STEP2_METRIC_COLUMNS))

    def test_step4_chromosome_accumulator_is_not_shadowed(self):
        tree = ast.parse((ROOT / 'scripts' / 'step4_summarize.py').read_text(),
                         filename='step4_summarize.py')
        stored_names = {node.id for node in ast.walk(tree)
                        if isinstance(node, ast.Name) and isinstance(node.ctx, ast.Store)}
        self.assertIn('chromosome_rows', stored_names)
        self.assertNotIn('chrom', stored_names)

if __name__ == '__main__':unittest.main()
