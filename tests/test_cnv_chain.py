"""Execute CNV boundary commands with synthetic files and mocked external tools."""
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile
import textwrap
import unittest

ROOT = Path(__file__).resolve().parents[1]
BASH = r'C:\Program Files\Git\bin\bash.exe' if os.name == 'nt' else shutil.which('bash')
HEADER = 'chromosome\tstart\tend\tgene\tlog2\tdepth\tweight\n'


def command(task, filename, inputs):
    source = (ROOT / filename).read_text(encoding='utf-8')
    block = re.search(r'task ' + re.escape(task) + r'\s*\{.*?command <<<(.*?)>>>', source, re.S).group(1)
    return textwrap.dedent(re.sub(r'~\{([^}]+)\}', lambda match: str(inputs[match[1]]), block))


class RegionBoundaryTest(unittest.TestCase):
    def run_region(self, contents, *options):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source, output = root / 'input.cnr', root / 'nested' / 'region.txt'
            source.write_text(contents)
            result = subprocess.run([sys.executable, str(ROOT / 'scripts/merge_bins_for_cnv.py'),
                                     '-i', str(source), '-o', str(output), *options], capture_output=True, text=True)
            return result, output.read_text() if output.exists() else ''

    def test_valid_zero_weight_target_and_antitarget_bins(self):
        source = HEADER + '1\t0\t100\tGENE\t-1\t30\t0\n1\t20000\t20100\tAntitarget\t1\t40\t0\n'
        result, output = self.run_region(source, '--bin-size', '20000', '--keep-antitarget')
        self.assertEqual(result.returncode, 0, result.stderr)
        rows = [row.split('\t') for row in output.splitlines()[1:]]
        self.assertEqual([row[3] for row in rows], ['DEL', 'DUP'])
        self.assertEqual([row[7] for row in rows], ['1.00', '4.00'])
        self.assertNotIn('nan', output.lower())

    def test_empty_nonfinite_malformed_and_invalid_options_fail_cleanly(self):
        normal = HEADER + '1\t0\t100\tGENE\t0\t30\t1\n'
        cases = [
            (HEADER, ['--bin-size', '20000']),
            (normal.replace('\t0\t30\t1', '\tnan\t30\t1'), ['--bin-size', '20000']),
            (normal.replace('\t0\t30\t1', '\tinf\t30\t1'), ['--bin-size', '20000']),
            (normal.replace('\t0\t30\t1', '\t2000\t30\t1'), ['--bin-size', '20000']),
            (normal.replace('\t30\t1', '\t30\t-1'), ['--bin-size', '20000']),
            ('bad header\n', ['--bin-size', '20000']),
            (HEADER + '1\t0\n', ['--bin-size', '20000']),
            (normal, ['--bin-size', '0']), (normal, ['--bin-count', '-1']),
            (normal, ['--bin-count', '1', '--dup-threshold', 'nan']),
            (normal, ['--bin-size', '20000', '--del-threshold', '3', '--dup-threshold', '2']),
        ]
        for contents, options in cases:
            with self.subTest(options=options, contents=contents):
                result, output = self.run_region(contents, *options)
                self.assertNotEqual(result.returncode, 0)
                self.assertIn('error:', result.stderr)
                self.assertNotIn('Traceback', result.stderr)
                self.assertEqual(output, '')

    def test_count_merging_and_default_antitarget_filter(self):
        source = HEADER + '1\t0\t100\tGENE\t0\t30\t1\n1\t100\t200\tAntitarget\t0\t30\t1\n'
        result, output = self.run_region(source, '--bin-count', '2')
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(len(output.splitlines()), 2)


class WDLBoundaryTest(unittest.TestCase):
    def shell(self, body, root, env=None):
        return subprocess.run([BASH, '-c', body], cwd=root, env=env, capture_output=True, text=True)

    def test_coverage_stops_at_first_failure_and_keeps_paths_with_spaces(self):
        for fail, expected in [(True, 17), (False, 0)]:
            with self.subTest(fail=fail), tempfile.TemporaryDirectory() as temporary:
                root = Path(temporary)
                body = command('CNVKitCoverage', 'tasks/cnvkit.wdl', {
                    'prefix': 'sample name', 'bam': 'input bam.bam', 'target_bed': 'target bed',
                    'antitarget_bed': 'anti bed', 'threads': 8,
                })
                fake = '''cnvkit.py() {
                    echo "$3" >> trace.txt
                    if [ "$3" = "target bed" ] && [ "$FAIL" = 1 ]; then return 17; fi
                    while [ "$#" -gt 0 ]; do if [ "$1" = -o ]; then shift; : > "$1"; break; fi; shift; done
                }
                '''
                result = self.shell(fake + body, root, dict(os.environ, FAIL='1' if fail else '0'))
                self.assertEqual(result.returncode, expected, result.stderr)
                self.assertEqual((root / 'trace.txt').read_text().splitlines(), ['target bed'] if fail else ['target bed', 'anti bed'])
                self.assertEqual((root / 'sample name.targetcoverage.cnn').exists(), not fail)

    def test_baseline_counts_are_validated_before_work(self):
        for left, right, valid in [(0, 0, False), (2, 1, False), (1, 1, True)]:
            with self.subTest(left=left, right=right), tempfile.TemporaryDirectory() as temporary:
                root = Path(temporary)
                body = command('ValidateBaselineInputs', 'tasks/germline.wdl', {'read_1_count': left, 'read_2_count': right})
                result = self.shell(body, root)
                self.assertEqual(result.returncode == 0, valid)
                self.assertEqual((root / 'sample_count.txt').exists(), valid)
        wdl = (ROOT / 'baseline.wdl').read_text(encoding='utf-8')
        self.assertIn('Int num_samples = ValidateInputs.sample_count', wdl)
        self.assertIn('scatter (i in range(num_samples))', wdl)

    def annotation(self, mode='ok', input_type='DEL'):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / 'input with spaces.txt'
            # A small gene-shaped table; the adapter must preserve every column.
            source.write_text('chromosome\tstart\tend\tcnv_call\tgene\n'
                              f'1\t0\t100\t{input_type}\tGENE1\n1\t100\t200\tDUP\tGENE2\n')
            fake = root / 'fake annotator.py'
            fake.write_text('''import argparse,os,sys
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('input');p.add_argument('-d');p.add_argument('-g');p.add_argument('-f');p.add_argument('-o');a=p.parse_args()
mode=os.environ['MODE']
if mode=='fail':sys.exit(23)
rows=Path(a.input).read_text().splitlines()[1:]
Path('seen.txt').write_text('\\n'.join(r.split('\\t')[3] for r in rows))
assert all(r.split('\\t')[3] in ['Deletion','Amplification','Normal'] for r in rows)
if mode=='drop':rows=rows[:-1]
if mode=='reverse':rows.reverse()
if mode=='mutate':rows[0]=rows[0].replace('GENE1','OTHER')
header='Chromosome\\tStart\\tEnd\\tCol4\\tCol5\\tISCN'
Path(a.o).write_text(header+'\\n'+'\\n'.join(r+'\\tannotation' for r in rows)+'\\n')
''')
            body = command('CNVAnno', 'tasks/cnvanno.wdl', {'cnv_bed': source.as_posix(), 'prefix': 'sample name', 'assembly': 'GRCh38'})
            body = body.replace('/app/CNVAnno/cnvanno.py', fake.as_posix())
            python = '"' + sys.executable.replace('\\', '/') + '"'
            body = body.replace('python - ', python + ' - ', 1)
            output = root / 'sample name.cnvanno.txt'
            output.write_text('previous result')
            result = self.shell(body, root, dict(os.environ, MODE=mode))
            return result, output.read_text(), (root / 'seen.txt').read_text() if (root / 'seen.txt').exists() else ''

    def test_annotation_normalizes_input_and_restores_original_types(self):
        result, output, seen = self.annotation()
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(seen.splitlines(), ['Deletion', 'Amplification'])
        self.assertEqual([row.split('\t')[3] for row in output.splitlines()[1:]], ['DEL', 'DUP'])
        self.assertIn('GENE1', output)

    def test_annotation_failure_or_changed_records_never_replace_valid_output(self):
        for mode in ['fail', 'drop', 'reverse', 'mutate']:
            with self.subTest(mode=mode):
                result, output, _ = self.annotation(mode)
                self.assertNotEqual(result.returncode, 0)
                self.assertEqual(output, 'previous result')
                if mode == 'fail':
                    self.assertEqual(result.returncode, 23)

    def test_sparse_cnr_to_gene_annotation_keeps_record_and_absolute_cn(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            cnr, gene = root / 'input.cnr', root / 'gene.txt'
            cnr.write_text(HEADER + ''.join(
                f'1\t{start}\t{start+100}\tGENE\t{value}\t30\t1\n'
                for start, value in zip([0, 100, 200, 100000, 100100, 100200], [-1.1, -1, -.9, -.1, 0, .1])))
            result = subprocess.run([sys.executable, str(ROOT / 'scripts/gene_level_segment.py'),
                                     '-i', str(cnr), '-o', str(gene)], capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            fake = root / 'annotator.py'
            fake.write_text("import sys\nfrom pathlib import Path\nrows=Path(sys.argv[1]).read_text().splitlines()[1:]\n"
                            "assert len(rows)==1 and rows[0].split('\\t')[3]=='Deletion'\n"
                            "header=['Chromosome','Start','End']+[f'Col{i}' for i in range(4,20)]+['ISCN']\n"
                            "Path(sys.argv[sys.argv.index('-o')+1]).write_text('\\t'.join(header)+'\\n'+rows[0]+'\\tannotation\\n')\n")
            body = command('CNVAnno', 'tasks/cnvanno.wdl', {'cnv_bed': gene.as_posix(), 'prefix': 'annotated', 'assembly': 'GRCh38'})
            body = body.replace('/app/CNVAnno/cnvanno.py', fake.as_posix()).replace('python - ',
                '\"' + sys.executable.replace('\\', '/') + '\" - ', 1)
            result = self.shell(body, root)
            self.assertEqual(result.returncode, 0, result.stderr)
            lines = (root / 'annotated.cnvanno.txt').read_text().splitlines()
            self.assertEqual(len(lines), 2)
            row = lines[1].split('\t')
            self.assertEqual(row[3], 'DEL')
            self.assertEqual(row[11], '1.41')
            self.assertEqual(row[14:16], ['6', '2'])

    def test_annotation_unknown_type_fails_before_launch(self):
        result, output, seen = self.annotation(input_type='UNKNOWN')
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual(seen, '')
        self.assertEqual(output, 'previous result')


if __name__ == '__main__':
    unittest.main()
