"""Run the actual WDL command against synthetic VCF/gVCF inputs in local Docker."""
import os
from pathlib import Path
import re
import subprocess
import textwrap
import unittest

ROOT = Path(__file__).resolve().parents[1]
HEADER = """##fileformat=VCFv4.2
##contig=<ID=17,length=1000>
##FILTER=<ID=LowQual,Description="Low quality">
##INFO=<ID=END,Number=1,Type=Integer,Description="End">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allele depth">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Likelihoods">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}
"""

def record(pos, alt, call, qual=".", end=None):
    return f"17\t{pos}\t.\tA\t{alt}\t{qual}\tPASS\t{'.' if end is None else 'END='+str(end)}\tGT:DP:GQ:AD:PL\t{call}\n"

@unittest.skipUnless(os.environ.get("SCHEMA_TEST_DOCKER") == "1", "set SCHEMA_TEST_DOCKER=1 for real bcftools integration")
class ProbandFamilyTest(unittest.TestCase):
    def run_family(self, proband, father, mother, expected, parent_order=("father", "mother"), expected_site="24.8\tPASS", max_seconds=30, slow=False, expected_exit=0):
        source = (ROOT / "tasks/glnexus.wdl").read_text(encoding="utf-8")
        block = re.search(r"task ProbandFamilyVCF\s*\{.*?command <<<(.*?)>>>", source, re.S).group(1)
        inputs = {
            'prefix': 'fixture', 'proband_vcf': 'proband.vcf.gz',
            'sep=" " parent_gvcfs': " ".join(f"{s}.vcf.gz" for s in parent_order),
            'sep=" " parent_gvcf_tbis': " ".join(f"{s}.vcf.gz.tbi" for s in parent_order),
            'sep=" " sample_names': "proband father mother",
            'sep="," sample_names': "proband,father,mother",
            'ref_dir': '.', 'fasta': 'ref.fa', 'threads': '1', 'max_seconds': str(max_seconds),
        }
        command = textwrap.dedent(re.sub(r"~\{([^}]+)\}", lambda m: inputs[m[1]], block))
        setup = "set -euo pipefail\ncd /tmp\n"
        setup += "printf '>17\\n' > ref.fa\nprintf '%1000s\\n' ' ' | tr ' ' A >> ref.fa\nprintf '17\\t1000\\t4\\t1000\\t1001\\n' > ref.fa.fai\n"
        for name, rows in [("proband", proband), ("father", father), ("mother", mother)]:
            setup += f"cat > {name}.vcf <<'VCF'\n" + HEADER.format(sample=name) + rows + "VCF\n"
            setup += f"bgzip {name}.vcf\nbcftools index -t {name}.vcf.gz\n"
        if slow:
            setup += 'mkdir -p /tmp/slowbin\ncat > /tmp/slowbin/bcftools <<\'SLOW\'\n#!/bin/bash\nsleep 3\nexec /usr/bin/bcftools "$@"\nSLOW\nchmod +x /tmp/slowbin/bcftools\nexport PATH=/tmp/slowbin:$PATH\n'
        finish = "\nprintf 'RESULT_START\\n'\nbcftools query -f '%POS\\t%REF\\t%ALT[\\t%GT:%DP:%GQ:%AD]\\n' fixture.proband-family.vcf.gz\n"
        finish += "printf 'SITE_START\\n'\nbcftools query -f '%QUAL\\t%FILTER\\n' fixture.proband-family.vcf.gz\n"
        result = subprocess.run(["docker", "run", "--rm", "--network", "none", "-i", "schemabio/glnexus:v1.4.1", "bash", "-s"],
                                input=(setup + command + finish).encode(), capture_output=True, timeout=90)
        result.stdout = result.stdout.decode()
        result.stderr = result.stderr.decode()
        self.assertEqual(result.returncode, expected_exit, result.stderr + result.stdout)
        if expected_exit:
            self.assertIn(f"failed exit_code={expected_exit}", result.stderr)
            return
        body, site = result.stdout.split("RESULT_START\n", 1)[1].split("SITE_START\n", 1)
        self.assertEqual(body.strip().splitlines(), expected)
        self.assertEqual(site.strip().splitlines(), [expected_site] * len(expected))

    def test_singleton_retained_with_uncovered_father_and_reference_mother(self):
        child = record(100, "C", "0/1:45:22:23,22:25,0,25", "24.8")
        father = record(1, "<*>", "0/0:30:45:30,.:0,45,90", end=90) + record(110, "<*>", "0/0:30:45:30,.:0,45,90", end=200)
        mother = record(1, "<*>", "0/0:30:45:30,.:0,45,90", end=200)
        self.run_family(child, father, mother,
                        ["100\tA\tC\t0/1:45:22:23,22\t./.:.:.:.\t0/0:30:45:30,."])

    def test_parental_variant_and_parent_only_sites(self):
        child = record(100, "C", "0/1:45:22:23,22:25,0,25", "24.8")
        father = record(100, "C", "0/1:40:30:20,20:30,0,30", "30") + record(150, "G", "0/1:40:30:20,20:30,0,30", "30")
        mother = record(1, "<*>", "0/0:30:45:30,.:0,45,90", end=200)
        self.run_family(child, father, mother,
                        ["100\tA\tC\t0/1:45:22:23,22\t0/1:40:30:20,20\t0/0:30:45:30,."], parent_order=("mother", "father"))

    def test_missing_site_quality_and_no_call_remain_unknown(self):
        child = record(100, "C", "0/1:45:22:23,22:25,0,25")
        father = record(1, "<*>", "./.:0:0:0,0:0,0,0", end=200)
        mother = record(1, "<*>", "0/0:30:45:30,0:0,45,90", end=200)
        self.run_family(child, father, mother,
                        ["100\tA\tC\t0/1:45:22:23,22\t./.:0:0:0,.\t0/0:30:45:30,."], expected_site=".\tPASS")

    def test_parent_pass_cannot_override_child_filter(self):
        child = record(100, "C", "0/1:45:22:23,22:25,0,25", "24.8").replace("\tPASS\t", "\tLowQual\t")
        father = record(100, "C", "0/1:40:30:20,20:30,0,30", "90")
        mother = record(1, "<*>", "0/0:30:45:30,0:0,45,90", end=200)
        self.run_family(child, father, mother,
                        ["100\tA\tC\t0/1:45:22:23,22\t0/1:40:30:20,20\t0/0:30:45:30,."], expected_site="24.8\tLowQual")

    def test_long_reference_block_supports_separated_proband_sites_once(self):
        child = record(100, "C", "0/1:45:22:23,22:25,0,25", "24.8") + record(150, "G", "0/1:45:22:23,22:25,0,25", "24.8")
        parent = record(1, "<*>", "0/0:30:45:30,0:0,45,90", end=200)
        self.run_family(child, parent, parent, [
            "100\tA\tC\t0/1:45:22:23,22\t0/0:30:45:30,.\t0/0:30:45:30,.",
            "150\tA\tG\t0/1:45:22:23,22\t0/0:30:45:30,.\t0/0:30:45:30,.",
        ])

    def test_no_proband_variants_returns_empty_family_report(self):
        father = record(100, "C", "0/1:40:30:20,20:30,0,30", "30")
        mother = record(1, "<*>", "0/0:30:45:30,0:0,45,90", end=200)
        self.run_family("", father, mother, [])

    def test_parent_indel_is_normalized_before_locus_selection(self):
        child = record(100, "A", "0/1:45:22:23,22:25,0,25", "24.8").replace("\tA\tA\t", "\tAA\tA\t")
        father = record(120, "A", "0/1:40:30:20,20:30,0,30", "30").replace("\tA\tA\t", "\tAA\tA\t")
        mother = record(1, "<*>", "0/0:30:45:30,0:0,45,90", end=200)
        self.run_family(child, father, mother,
                        ["1\tAA\tA\t0/1:45:22:23,22\t0/1:40:30:20,20\t0/0:30:45:30,."])

    def test_whole_task_timeout_is_explicit(self):
        child = record(100, "C", "0/1:45:22:23,22:25,0,25", "24.8")
        parent = record(1, "<*>", "0/0:30:45:30,0:0,45,90", end=200)
        self.run_family(child, parent, parent, [], max_seconds=1, slow=True, expected_exit=124)

    def test_multiallelic_child_and_parent_with_different_alt(self):
        child = record(100, "C,G", "1/2:45:22:0,20,25:90,60,80,60,0,80", "24.8")
        father = record(100, "C,<*>", "0/1:40:30:20,20,0:30,0,30,90,90,90", "30")
        mother = record(100, "G,<*>", "0/1:40:30:20,20,0:30,0,30,90,90,90", "30")
        self.run_family(child, father, mother, [
            "100\tA\tC\t1/0:45:22:0,20\t0/1:40:30:20,20\t0/0:40:30:20,.",
            "100\tA\tG\t0/1:45:22:0,25\t0/0:40:30:20,.\t0/1:40:30:20,20",
        ])

if __name__ == "__main__":
    unittest.main()
