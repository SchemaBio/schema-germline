"""CBS boundary and CLI regressions using synthetic data only."""
import importlib.util
import math
from pathlib import Path
import random
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / 'scripts' / 'gene_level_segment.py'
spec = importlib.util.spec_from_file_location('gene_level_segment', SCRIPT)
cnv = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = cnv
spec.loader.exec_module(cnv)
HEADER = 'chromosome\tstart\tend\tgene\tlog2\tdepth\tweight\n'


def bins(values, weight=1.0):
    return [cnv.Bin('1', i * 100, (i + 1) * 100, 'GENE|TX|ENSG|1|+|p', v, 30.0, weight)
            for i, v in enumerate(values)]


class CBSBoundaryTest(unittest.TestCase):
    def segment(self, values, min_bins=3, weights=None, **kwargs):
        return cnv.cbs_segment_recursive(values, list(range(len(values))),
                                         [1.0] * len(values) if weights is None else weights,
                                         min_bins=min_bins, **kwargs)

    def test_incident_singleton_right_variance_is_not_evaluated(self):
        self.assertEqual(cnv.calculate_cbs_statistic([0, 0, 0, 1], 2), 0)
        with patch.object(cnv, 'calculate_cbs_statistic', side_effect=AssertionError('invalid split evaluated')):
            self.assertEqual(self.segment([0, 0, 0, 0, 0, 10]), [(0, 5)])

    def test_statistic_rejects_empty_singleton_and_out_of_range_splits(self):
        for values in ([], [1], [1, 2], [1, 2, 3], [1, 2, 3, 4]):
            for k in (-2, -1, 0, len(values) - 2, len(values) - 1, len(values), 100):
                with self.subTest(values=values, k=k):
                    self.assertEqual(cnv.calculate_cbs_statistic(values, k), 0)

    def test_valid_statistic_preserves_existing_formula(self):
        self.assertAlmostEqual(cnv.calculate_cbs_statistic([0, 2, 4, 6], 1), 2.0)

    def test_exact_minimum_bins_on_both_sides_can_split(self):
        self.assertEqual(self.segment([-1.1, -1.0, -0.9, 0.9, 1.0, 1.1]), [(0, 2), (3, 5)])

    def test_empty_singleton_short_and_constant_inputs(self):
        self.assertEqual(self.segment([]), [])
        self.assertEqual(cnv.run_cbs_on_chromosome([]), [])
        for n in range(1, 12):
            self.assertEqual(self.segment([0.0] * n), [(0, n - 1)])
        segment = cnv.run_cbs_on_chromosome(bins([0.2]))[0]
        self.assertEqual(segment.log2_std, 0)
        self.assertEqual(segment.bin_count, 1)

    def test_zero_weights_keep_existing_unweighted_fallback(self):
        values = [-1.1, -1.0, -0.9, 0.9, 1.0, 1.1]
        self.assertEqual(self.segment(values, weights=[0] * 6), self.segment(values))
        self.assertTrue(all(math.isfinite(s.log2_mean) for s in cnv.run_cbs_on_chromosome(bins(values, 0))))

    def test_maximum_depth_stops_recursion(self):
        self.assertEqual(self.segment([-1] * 6 + [1] * 6, max_depth=0), [(0, 11)])

    def test_invalid_parameters_lengths_and_values_are_explicit_errors(self):
        for alpha in (0, -1, 2, math.nan, math.inf):
            with self.assertRaisesRegex(ValueError, 'alpha'):
                self.segment([0] * 6, alpha=alpha)
        for minimum in (0, 1, -1, 2.5):
            with self.assertRaisesRegex(ValueError, 'min_bins'):
                self.segment([0] * 6, min_bins=minimum)
        with self.assertRaisesRegex(ValueError, 'equal lengths'):
            self.segment([0] * 6, weights=[1])
        for value in (math.nan, math.inf, -math.inf):
            with self.assertRaisesRegex(ValueError, 'finite'):
                self.segment([value])
            with self.assertRaisesRegex(ValueError, 'finite'):
                cnv.run_cbs_on_chromosome(bins([value]))
        with self.assertRaisesRegex(ValueError, 'weight'):
            self.segment([0], weights=[-1])

    def test_synthetic_partitions_cover_all_bins_without_gaps_or_overlaps(self):
        rng = random.Random(20261031)
        for n in range(1, 100):
            values = [rng.gauss(-1 if i < n // 2 else 1, 0.2) for i in range(n)]
            minimum = 2 + n % 5
            segments = self.segment(values, min_bins=minimum)
            self.assertEqual([i for start, end in segments for i in range(start, end + 1)], list(range(n)))
            if len(segments) > 1:
                self.assertTrue(all(end - start + 1 >= minimum for start, end in segments))


    def test_sparse_gene_exons_are_not_dropped_after_valid_segmentation(self):
        values = [-1.1, -1.0, -0.9, 0.9, 1.0, 1.1]
        records = bins(values)
        for record, start in zip(records, [0, 100, 200, 100000, 100100, 100200]):
            record.start, record.end = start, start + 100
        segments = cnv.run_cbs_segmentation(records)
        self.assertEqual(len(segments['1']), 2)
        genes = cnv.gene_level_cnv(records, segments)
        self.assertEqual(len(genes), 1)
        self.assertEqual(genes[0].bin_count, 6)
        self.assertEqual(genes[0].segment_count, 2)
        self.assertLess(genes[0].coverage_ratio, 0.1)

    def test_missing_segment_coverage_is_an_explicit_contract_error(self):
        with self.assertRaisesRegex(ValueError, 'cover all bins'):
            cnv.gene_level_cnv(bins([0, 0, 0]), {'1': []})


class CNRInputTest(unittest.TestCase):
    def read(self, contents):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / 'input.cnr'
            path.write_text(contents)
            return cnv.read_cnr_file(str(path))

    def test_blank_lines_and_antitargets_are_allowed(self):
        contents = HEADER + '\n1\t0\t100\tAntitarget\t0\t30\t1\n1\t100\t200\tGENE\t0\t30\t0\n'
        self.assertEqual(len(self.read(contents)), 1)

    def test_invalid_rows_include_source_line_and_reason(self):
        for row in ('1\t0', '1\t0\t100\tGENE\tnan\t30\t1',
                    '1\t0\t100\tGENE\t0\tinf\t1', '1\t0\t100\tGENE\t0\t30\t-1',
                    '1\t100\t0\tGENE\t0\t30\t1', '1\tX\t100\tGENE\t0\t30\t1'):
            with self.subTest(row=row), self.assertRaisesRegex(ValueError, r'input.cnr:2:'):
                self.read(HEADER + row + '\n')
        with self.assertRaisesRegex(ValueError, ':1: invalid'):
            self.read('wrong header\n')

    def cli(self, contents, *options):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / 'input.cnr'
            source.write_text(contents)
            output = root / 'genes.tsv'
            segments = root / 'segments.tsv'
            result = subprocess.run([sys.executable, str(SCRIPT), '-i', str(source), '-o', str(output),
                                     '--segments', str(segments), *options], capture_output=True, text=True)
            return result, output.read_text() if output.exists() else '', segments.read_text() if segments.exists() else ''

    def test_cli_incident_shape_produces_finite_gene_and_segment_outputs(self):
        values = [0] * 5 + [10]
        contents = HEADER + ''.join(f'1\t{i * 100}\t{(i+1) * 100}\tGENE\t{v}\t30\t1\n' for i, v in enumerate(values))
        result, output, segments = self.cli(contents)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(len(output.splitlines()), 2)
        self.assertEqual(len(segments.splitlines()), 2)
        self.assertNotIn('nan', output.lower())
        self.assertNotIn('inf', output.lower())

    def test_cli_empty_or_invalid_input_fails_without_publishing_results(self):
        for contents, options in ((HEADER, []), (HEADER + '1\t0\t100\tGENE\tnan\t30\t1\n', []),
                                  (HEADER + '1\t0\t100\tGENE\t2000\t30\t1\n', []),
                                  (HEADER, ['--cn-del-threshold', 'nan']),
                                  (HEADER, ['--cn-del-threshold', '3', '--cn-dup-threshold', '2']),
                                  (HEADER, ['--min-bins', '1']), (HEADER, ['--alpha', 'nan'])):
            result, output, segments = self.cli(contents, *options)
            self.assertNotEqual(result.returncode, 0)
            self.assertIn('error:', result.stderr)
            self.assertNotIn('Traceback', result.stderr)
            self.assertEqual((output, segments), ('', ''))


if __name__ == '__main__':
    unittest.main()
