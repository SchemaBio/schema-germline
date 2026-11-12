"""Real GATK regressions for normalization reordering and partial output."""
import json
import os
from pathlib import Path
import re
import subprocess
import textwrap
import unittest

ROOT = Path(__file__).resolve().parents[1]

SETUP = r"""set -euo pipefail
cd /tmp
python3 - <<'FIXTURE'
from pathlib import Path
chroms=['1','2','3','4','17']
with open('ref.fa','w') as fasta, open('ref.fa.fai','w') as fai:
    for chrom in chroms:
        fasta.write('>'+chrom+'\n'); offset=fasta.tell()
        sequence=list('A'*1000)
        if chrom=='3': sequence[101]='C'
        fasta.write(''.join(sequence)+'\n')
        fai.write(f'{chrom}\t1000\t{offset}\t1000\t1001\n')
Path('ref.dict').write_text('@HD\tVN:1.6\tSO:coordinate\n'+''.join(f'@SQ\tSN:{c}\tLN:1000\n' for c in chroms))
header='##fileformat=VCFv4.2\n'+''.join(f'##contig=<ID={c},length=1000>\n' for c in chroms)+'##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tfather\tmother\tproband\n'
records=[('1',10,'A','C,G','1/2'),('2',10,'A','C','0/1'),('3',50,'A','C','0/1'),('3',101,'AC','ACT,A','1/2'),('4',10,'A','C','0/1'),('17',10,'A','C','0/1')]
Path('input.vcf').write_text(header+''.join(f'{c}\t{pos}\t.\t{ref}\t{alt}\t50\tPASS\t.\tGT\t0/0\t0/0\t{gt}\n' for c,pos,ref,alt,gt in records))
FIXTURE
"""
FINISH = r"""
code=$?
export TASK_EXIT=$code
python3 - <<'RESULT'
import gzip,json,os
from pathlib import Path
from collections import Counter
p=Path('fixture.left.vcf.gz'); rows=[]
if p.exists():
    with gzip.open(p,'rt') as f: rows=[line.rstrip().split('\t') for line in f if not line.startswith('#')]
print('RESULT '+json.dumps({'exit':int(os.environ['TASK_EXIT']),'published':p.exists(),'rows':len(rows),'chromosomes':dict(Counter(row[0] for row in rows)),'sites':[[r[0],int(r[1]),r[3],r[4]] for r in rows],'index':Path(str(p)+'.tbi').exists()}))
RESULT
exit "$code"
"""

@unittest.skipUnless(os.environ.get("SCHEMA_TEST_DOCKER") == "1", "set SCHEMA_TEST_DOCKER=1 for GATK integration")
class LeftAlignOutputTest(unittest.TestCase):
    def run_task(self, command=None, failure=None):
        if command is None:
            source = (ROOT / "tasks/gatk.wdl").read_text(encoding="utf-8")
            command = re.search(r"task LeftAlignAndTrimVariants\s*\{.*?command <<<(.*?)>>>",source,re.S).group(1)
        values={"ref_dir":".","fasta":"ref.fa","vcf":"input.vcf","prefix":"fixture"}
        command=textwrap.dedent(re.sub(r"~\{([^}]+)\}",lambda m:values[m[1]],command))
        setup=SETUP
        if failure:
            wrapper = '#!/bin/bash\nif [ "$1" = "'+failure+'" ]; then echo "INJECTED_TOOL_FAILURE" >&2; exit 42; fi\nexec "$REAL_GATK" "$@"\n'
            setup += "export REAL_GATK=$(command -v gatk)\nmkdir -p /tmp/mockbin\ncat > /tmp/mockbin/gatk <<'WRAPPER'\n"+wrapper+"WRAPPER\nchmod +x /tmp/mockbin/gatk\nexport PATH=/tmp/mockbin:$PATH\n"
        # Model MiniWDL's ordinary shell. The task must enforce its own errors.
        script=setup+"set +e\n(\n"+command+"\n)\n"+FINISH
        result=subprocess.run(["docker","run","--rm","--network","none","-i","--cpus","2","--memory","3g","schemabio/gatk:4.6.2.0","bash","-s"],
                              input=script.encode(),capture_output=True,timeout=90)
        lines=[line[7:] for line in result.stdout.decode().splitlines() if line.startswith("RESULT ")]
        self.assertEqual(len(lines),1,result.stderr.decode()[-5000:])
        return result,json.loads(lines[0])

    def test_multiallelic_trimming_is_sorted_and_all_chromosomes_survive(self):
        result,summary=self.run_task()
        self.assertEqual(result.returncode,0,result.stderr.decode()[-5000:])
        self.assertEqual(summary["chromosomes"],{"1":2,"2":1,"3":3,"4":1,"17":1})
        self.assertEqual(summary["rows"],8)
        self.assertTrue(summary["index"])
        self.assertIn(["3",101,"AC","A"],summary["sites"])
        self.assertIn(["3",102,"C","CT"],summary["sites"])
        chr3=[r[1] for r in summary["sites"] if r[0]=="3"]
        self.assertEqual(chr3,sorted(chr3))

    def test_normalization_failure_cannot_be_masked_by_indexing(self):
        result,summary=self.run_task(failure="LeftAlignAndTrimVariants")
        self.assertEqual(result.returncode,42)
        self.assertFalse(summary["published"])
        self.assertFalse(summary["index"])

    def test_sort_failure_does_not_publish_output(self):
        result,summary=self.run_task(failure="SortVcf")
        self.assertEqual(result.returncode,42)
        self.assertFalse(summary["published"])
        self.assertFalse(summary["index"])

if __name__ == "__main__":
    unittest.main()

