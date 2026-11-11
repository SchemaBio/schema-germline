"""Offline Docker benchmark using synthetic records only; --output is required."""
import argparse,json,re,subprocess,textwrap,time
from pathlib import Path

p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);p.add_argument('--sites',type=int,default=62500);p.add_argument('--block',type=int,default=100);p.add_argument('--run',action='store_true');p.add_argument('--name',default='baseline');p.add_argument('--staggered',action='store_true');a=p.parse_args()
if a.sites<0 or a.block<1 or not re.fullmatch(r'[A-Za-z0-9_-]+',a.name): p.error('sites must be nonnegative, block positive, and name alphanumeric')
root=a.output.resolve();root.mkdir(parents=True,exist_ok=True)
length=a.sites*1000+1000
header='''##fileformat=VCFv4.2
##contig=<ID=17,length=LENGTH>
##INFO=<ID=END,Number=1,Type=Integer,Description="End">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allele depth">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Likelihoods">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE
'''.replace('LENGTH',str(length))
if not (root/'fixture.json').exists():
 with (root/'ref.fa').open('w',newline='\n') as f:
  f.write('>17\n')
  for i in range(0,length,100): f.write(('ACGT'*25)[:min(100,length-i)]+'\n')
 (root/'ref.fa.fai').write_text(f'17\t{length}\t4\t100\t101\n')
 with (root/'proband.vcf').open('w',newline='\n') as f:
  f.write(header.replace('SAMPLE','proband'))
  for i in range(a.sites): f.write(f'17\t{i*1000+101}\t.\tA\tC\t24.8\tPASS\t.\tGT:DP:GQ:AD:PL\t0/1:45:22:23,22:25,0,25\n')
 for name in ['father','mother']:
  step=a.block+int(a.staggered and name=='mother')
  with (root/(name+'.vcf')).open('w',newline='\n') as f:
   f.write(header.replace('SAMPLE',name))
   for start in range(1,length,step):
    ref='ACGT'[(start-1)%4]
    f.write(f'17\t{start}\t.\t{ref}\t<*>\t.\tPASS\tEND={min(start+step-1,length)}\tGT:DP:GQ:AD:PL\t0/0:30:45:30,0:0,45,90\n')
 (root/'fixture.json').write_text(json.dumps({'proband_sites':a.sites,'parent_refblocks_each':(length+a.block-1)//a.block,'reference_bases':length}))
source=(Path(__file__).resolve().parents[1]/'tasks/glnexus.wdl').read_text(encoding='utf-8')
block=re.search(r'task ProbandFamilyVCF\s*\{.*?command <<<(.*?)>>>',source,re.S).group(1)
inputs={'prefix':'fixture','proband_vcf':'/bench/proband.vcf.gz','sep=" " parent_gvcfs':'/bench/father.vcf.gz /bench/mother.vcf.gz','sep=" " parent_gvcf_tbis':'/bench/father.vcf.gz.tbi /bench/mother.vcf.gz.tbi','sep=" " sample_names':'proband father mother','sep="," sample_names':'proband,father,mother','ref_dir':'/bench','fasta':'ref.fa','threads':'4','max_seconds':'120'}
command=textwrap.dedent(re.sub(r'~\{([^}]+)\}',lambda m:inputs[m[1]],block))
setup='set -euo pipefail\nfor s in proband father mother; do if [ ! -s /bench/$s.vcf.gz.tbi ]; then bgzip -c /bench/$s.vcf > /bench/$s.vcf.gz; bcftools index -t /bench/$s.vcf.gz; fi; done\nmkdir -p /bench/'+a.name+'\ncd /bench/'+a.name+'\nPS4=\'+${SECONDS}s \'\nset -x\n'
(root/(a.name+'.sh')).write_text(setup+command,encoding='utf-8',newline='\n')
print((root/'fixture.json').read_text(),flush=True)
if a.run:
 started=time.monotonic()
 with (root/(a.name+'.stdout')).open('wb') as out,(root/(a.name+'.stderr')).open('wb') as err:
  r=subprocess.run(['docker','run','--rm','--network','none','--cpus','4','--memory','8g','--mount',f'type=bind,source={root.resolve()},target=/bench','schemabio/glnexus:v1.4.1','bash','/bench/'+a.name+'.sh'],stdout=out,stderr=err,timeout=180)
 print(json.dumps({'run':a.name,'returncode':r.returncode,'elapsed_seconds':round(time.monotonic()-started,3)}),flush=True)
