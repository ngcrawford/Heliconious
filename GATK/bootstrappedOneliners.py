import shlex
import random
from subprocess import Popen, PIPE
def get_subset_vcf(chrm, start, stop):
    base_dir = "/Users/MullenLab/Desktop/Grad_Students/Nick/butterfly_practice"
    cli = """java -Xmx6g -jar /Users/MullenLab/Source/gatk/dist/GenomeAnalysisTK.jar \
      -R {3}/Butterfly.merge.scaffolds.fa  \
      -T SelectVariants \
      --variant {3}/32_butterflies_vcfs/32.butterflies.good_BGI_snps.combined.vcf \
      -L {0}:{1}-{2}""".format(chrm, start, stop, base_dir)

    cli_parts = shlex.split(cli)
    vcf = Popen(cli_parts, stdin=PIPE, stderr=PIPE, stdout=PIPE).communicate()[0]
    return vcf

def generate_bootstraps(chrm, chrm_len, window_size, numb_of_reps):

    start_sites = [random.randint(0, chrm_len-1000) for item in range(numb_of_reps)]
    replicates = []
    for start in start_sites:
        vcf = get_subset_vcf(chrm, start, start+window_size)
        print chrm, start, start+window_size
        vcf = vcf.split("\n")
        snps = [line.split("\t") for line in vcf if line.startswith(chrm)]
        replicates.append(snps)

    return replicates

chrm = 'Yb_superscaffold'
chrm_len = 1333113
window_size = 10000
numb_of_reps = 3

replicates = generate_bootstraps(chrm, chrm_len, window_size, numb_of_reps)