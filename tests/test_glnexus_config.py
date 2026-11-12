"""Exercise the GLNexus WDL command with its actual v1.4.1 binary."""
import os
from pathlib import Path
import re
import subprocess
import textwrap
import unittest

ROOT = Path(__file__).resolve().parents[1]
HEADER = """##fileformat=VCFv4.2
##contig=<ID=17,length=1000>
##INFO=<ID=END,Number=1,Type=Integer,Description="End">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">
##FORMAT=<ID=MIN_DP,Number=1,Type=Integer,Description="Minimum reference depth">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Likelihoods">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}
"""


class TrioRouteTest(unittest.TestCase):
    def test_joint_calling_is_active_and_proband_task_is_retained(self):
        workflow = (ROOT / "trio.wdl").read_text(encoding="utf-8")
        tasks = (ROOT / "tasks/glnexus.wdl").read_text(encoding="utf-8")
        self.assertNotIn("call GLNEXUS.ProbandFamilyVCF", workflow)
        self.assertIn("task ProbandFamilyVCF {", tasks)
        split = re.search(r"call GERMLINE.SplitVcf as SplitVcfHap\s*\{(.*?)\}", workflow, re.S).group(1)
        self.assertIn("vcf = GLNexus.out_vcf", split)


@unittest.skipUnless(os.environ.get("SCHEMA_TEST_DOCKER") == "1", "set SCHEMA_TEST_DOCKER=1 for GLNexus integration")
class GLNexusAQ20Test(unittest.TestCase):
    def run_case(self, aq=24, covered=True):
        source = (ROOT / "tasks/glnexus.wdl").read_text(encoding="utf-8")
        block = re.search(r"task GLNexus\s*\{.*?command <<<(.*?)>>>", source, re.S).group(1)
        names = ["proband", "father", "mother"]
        inputs = {
            "prefix": "fixture", "threads": "1",
            'sep=" " gvcf': " ".join(f"/tmp/{name}.g.vcf.gz" for name in names),
            'sep=" " gvcf_tbi': " ".join(f"/tmp/{name}.g.vcf.gz.tbi" for name in names),
        }
        command = textwrap.dedent(re.sub(r"~\{([^}]+)\}", lambda m: inputs[m[1]], block))
        setup = "set -euo pipefail\ncd /tmp\n"
        for name in names:
            if name == "proband":
                row = f"17\t100\t.\tA\tC,<*>\t24.8\tPASS\t.\tGT:DP:GQ:AD:PL\t0/1:45:{aq}:23,22,0:{aq},0,80,90,90,90\n"
            else:
                pos = 1 if covered else 500
                row = f"17\t{pos}\t.\tA\t<*>\t.\tPASS\tEND=1000\tGT:DP:MIN_DP:GQ:AD:PL\t0/0:30:30:60:30,0:0,60,90\n"
            setup += f"cat > {name}.g.vcf <<'VCF'\n" + HEADER.format(sample=name) + row + "VCF\n"
            setup += f"bgzip {name}.g.vcf\nbcftools index -t {name}.g.vcf.gz\n"
        setup += "glnexus_cli --config DeepVariantWES --dir baseline.DB --threads 1 --mem-gbytes 1 proband.g.vcf.gz father.g.vcf.gz mother.g.vcf.gz > baseline.bcf\n"
        finish = """
printf 'BASE_CONFIG\\n'
bcftools view -h baseline.bcf | grep '^##GLnexusConfig='
printf 'BASE_ROWS\\n'
bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT[\\t%SAMPLE=%GT:GQ=%GQ:DP=%DP]\\n' baseline.bcf
printf 'CURRENT_CONFIG\\n'
bcftools view -h fixture.glnexus.vcf.gz | grep '^##GLnexusConfig='
printf 'CURRENT_ROWS\\n'
bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT[\\t%SAMPLE=%GT:GQ=%GQ:DP=%DP]\\n' fixture.glnexus.vcf.gz
"""
        result = subprocess.run(
            ["docker", "run", "--rm", "--network", "none", "-i", "--cpus", "1", "--memory", "2g", "schemabio/glnexus:v1.4.1", "bash", "-s"],
            input=(setup + command + finish).encode(), capture_output=True, timeout=60,
        )
        self.assertEqual(result.returncode, 0, result.stderr.decode()[-4000:])
        sections = {key: [] for key in ["BASE_CONFIG", "BASE_ROWS", "CURRENT_CONFIG", "CURRENT_ROWS"]}
        current = None
        for line in result.stdout.decode().splitlines():
            if line in sections:
                current = line
            elif current and line:
                sections[current].append(line)
        # Compare the binary's complete effective config, including defaults and
        # field liftover definitions. Exactly one value may differ.
        self.assertEqual(len(sections["BASE_CONFIG"]), 1)
        self.assertEqual(sections["CURRENT_CONFIG"], [
            sections["BASE_CONFIG"][0].replace("min_AQ1: 35", "min_AQ1: 20", 1)
        ])
        return sections

    def test_singleton_aq24_is_retained_with_parent_reference_calls(self):
        result = self.run_case()
        self.assertEqual(result["BASE_ROWS"], [])
        self.assertEqual(result["CURRENT_ROWS"], [
            "17\t100\tA\tC\tfather=0/0:GQ=60:DP=30\tmother=0/0:GQ=60:DP=30\tproband=0/1:GQ=24:DP=45"
        ])

    def test_uncovered_parents_remain_missing(self):
        result = self.run_case(covered=False)
        self.assertEqual(result["BASE_ROWS"], [])
        self.assertEqual(result["CURRENT_ROWS"], [
            "17\t100\tA\tC\tfather=./.:GQ=.:DP=.\tmother=./.:GQ=.:DP=.\tproband=0/1:GQ=24:DP=45"
        ])

    def test_aq19_is_not_promoted_to_unfiltered_output(self):
        self.assertEqual(self.run_case(aq=19)["CURRENT_ROWS"], [])


if __name__ == "__main__":
    unittest.main()
