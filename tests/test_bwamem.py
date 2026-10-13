"""Execute the WDL alignment command with mocked aligners and reference indexes."""
import os
from pathlib import Path
import re
import shutil
import subprocess
import tempfile
import unittest

BASH = shutil.which('bash')
if os.name == 'nt':
    BASH = r'C:\Program Files\Git\bin\bash.exe'


class BwaAlignCommandTest(unittest.TestCase):
    def run_alignment(self, indexes, *, threads=8, avx=True, memory=33554432, fail=False, minibwa=True):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            ref = root / 'ref'
            ref.mkdir()
            for suffix in indexes:
                (ref / ('reference.fa.' + suffix)).write_text('index')
            tools = root / 'bin'
            tools.mkdir()
            for mapper in ['minibwa', 'bwa']:
                if mapper == 'minibwa' and not minibwa:
                    continue
                tool = tools / mapper
                tool.write_text(f'#!/bin/bash\necho {mapper} >> "$TRACE"\n' + ('exit 3\n' if fail else 'echo SAM\n'))
                tool.chmod(0o755)
            samtools = tools / 'samtools'
            samtools.write_text('#!/bin/bash\nif [ "$1" = sort ]; then cat > result.sorted.bam; else touch result.sorted.bam.bai; fi\n')
            samtools.chmod(0o755)
            wdl = (Path(__file__).parents[1] / 'tasks/bwamem.wdl').read_text(encoding='utf-8')
            command = wdl.split('command <<<', 1)[1].split('>>>', 1)[0]
            substitutions = {'actual_threads': str(threads), 'ref_dir': ref.as_posix(), 'ref_fasta_name': 'reference.fa', 'read_1': 'read1.fq', 'read_2': 'read2.fq', 'prefix': 'result'}
            command = re.sub(r'~\{([^}]+)\}', lambda match: substitutions[match[1]], command)
            preamble = f'''grep() {{
    case "$*" in
        *cpuinfo*) return {0 if avx else 1} ;;
        *meminfo*) echo 'MemTotal: {memory} kB' ;;
        *) command grep "$@" ;;
    esac
}}
'''
            env = dict(os.environ, PATH=tools.as_posix() + os.pathsep + os.environ['PATH'], TRACE=(root / 'trace').as_posix())
            if os.name == 'nt':
                env['PATH'] = tools.as_posix() + os.pathsep + r'C:\Program Files\Git\usr\bin' + os.pathsep + os.environ['PATH']
            result = subprocess.run([BASH, '-c', preamble + command], cwd=root, env=env, capture_output=True, text=True)
            trace = (root / 'trace').read_text().strip() if (root / 'trace').exists() else ''
            return result, trace

    def test_minibwa_is_preferred_even_with_eight_threads(self):
        result, mapper = self.run_alignment(['l2b', 'mbw', 'amb', 'ann', 'bwt', 'pac', 'sa'])
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(mapper, 'minibwa')

    def test_standard_bwa_fallback_requires_complete_indexes(self):
        result, mapper = self.run_alignment(['amb', 'ann', 'bwt', 'pac', 'sa'])
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(mapper, 'bwa')

    def test_incompatible_cpu_or_memory_can_use_complete_bwa(self):
        for options in [{'avx': False}, {'memory': 8388608}, {'minibwa': False}]:
            with self.subTest(options=options):
                result, mapper = self.run_alignment(['l2b', 'mbw', 'amb', 'ann', 'bwt', 'pac', 'sa'], **options)
                self.assertEqual(result.returncode, 0, result.stderr)
                self.assertEqual(mapper, 'bwa')

    def test_incomplete_indexes_fail_before_alignment(self):
        result, mapper = self.run_alignment(['l2b', 'amb', 'ann', 'bwt', 'pac'])
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual(mapper, '')
        self.assertIn('missing or empty', result.stderr)

    def test_failed_alignment_is_not_hidden_by_samtools(self):
        result, mapper = self.run_alignment(['l2b', 'mbw'], fail=True)
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual(mapper, 'minibwa')


if __name__ == '__main__':
    unittest.main()
