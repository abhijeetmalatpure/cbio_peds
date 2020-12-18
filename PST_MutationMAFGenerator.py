import os
from os.path import expanduser, join, sep
import subprocess
from configparser import ConfigParser
import time, datetime
import logging
import socket

logger = None

# Create log directory
if not os.path.isdir('log'):
    os.mkdir('log')

# Before we do anything, set up logging.
epoch_time = time.mktime(datetime.datetime.now().timetuple())
log_filename = os.path.join(f"log/mafconverter_{str(int(epoch_time))}.log")

print(f'Writing to {log_filename}')
logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s %(message)s", filename=log_filename)
logger = logging.getLogger()

logger.info(f"Starting download on {socket.gethostname()} with process ID: {os.getpid()}")
# Each list within vcf_files has-
# 1. File path
# 2. Sample name to be used - "filename"
# 3. Sample name found in vcf file
# The output maf and enhanced maf will have names
# filename.somatic.enhanced.vcf, filename.somatic.maf and filename.somatic.enhanced.maf (or germline)
vcfs = []
header = []

for root, dirs, files in os.walk("/N/slate/abhmalat/peds_pst/ashion/files/subjects"):
        for file in files:
            if file.endswith(".somatic.vcf") and root.split(sep)[-1] == 'normalized':
                with open(join(root, file), 'r') as content:
                    for line in content:
                        line = line.strip()
                        if line.startswith("#CHROM\t"):
                            vcfs.append([join(root, file),
                                         file.split('.')[0],
                                         line.split("\t")[-1]])
                            break

logger.info(f'Somatic files found:{len(vcfs)} ')
logger.info(vcfs)

[os.makedirs(join(expanduser("~"), 'tmp', (vcf[1] + 'somatic')),
             exist_ok=True) for vcf in vcfs]

# Load samtools, bcftools
config = ConfigParser()
config.read('download.cfg')
libraries = dict(config.items('libraries'))

os.environ['LOADEDMODULES'] = libraries['loadedmodules']
os.environ['PATH'] = libraries['path']
os.environ['_LMFILES_'] = libraries['lmfiles']
os.environ['LD_LIBRARY_PATH'] = libraries['ld_library_path']

toolspath = join(expanduser("~"), "cbio_tools", "vcf2maf")
vcfpath = "/N/slate/abhmalat/peds_pst/ashion/vcf2mafConversion/vcf"
mafpath = "/N/slate/abhmalat/peds_pst/ashion/vcf2mafConversion/maf"
enhancedpath = "/N/slate/abhmalat/peds_pst/ashion/vcf2mafConversion/enhanced"

subprocess.Popen(['ls'], stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()

success = 0
for vcf in vcfs:
    vcf2vcf = [
        'perl', join(toolspath, 'vcf2vcf.pl'),
        '--input-vcf', vcf[0],
        '--output-vcf', join(vcfpath, (vcf[1] + '.somatic.enhanced.vcf')),
        '--vcf-tumor-id', vcf[2],
        '--new-tumor-id', vcf[1],
        '--ref-fasta',
        join(expanduser("~"),
             '.vep/homo_sapiens/102_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz')]

    vcf2maf = [
            'perl', join(toolspath, 'vcf2maf.pl'),
            '--input-vcf', join(vcfpath, (vcf[1] + '.somatic.enhanced.vcf')),
            '--output-maf', join(mafpath, (vcf[1] + '.somatic.maf')),
            '--vep-path', join(toolspath, 'vep_hg19'),
            '--vep-data', join(expanduser("~"), '.vep'),
            '--vep-forks', '5',
            '--tmp-dir', join(expanduser("~"), 'tmp', (vcf[1] + 'somatic')),
            '--vcf-tumor-id', vcf[1],
            '--tumor-id', vcf[1],
            '--ncbi-build', 'GRCh37',
            '--ref-fasta',
            join(expanduser("~"),
                 '.vep/homo_sapiens/102_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz')
        ]

    maf2maf = [
                'perl', join(toolspath, 'maf2maf.pl'),
                '--input-maf', join(mafpath, (vcf[1] + '.somatic.maf')),
                '--output-maf', join(enhancedpath, (vcf[1] + '.somatic.enhanced.maf')),
                '--vep-path', join(toolspath, 'vep_hg19'),
                '--vep-data', join(expanduser("~"), '.vep'),
                '--vep-forks', '5',
                '--tmp-dir', join(expanduser("~"), 'tmp', (vcf[1] + 'somatic')),
                '--ncbi-build', 'GRCh37',
                '--ref-fasta',
                join(expanduser("~"),
                     '.vep/homo_sapiens/102_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz')
            ]


    subprocess.call(vcf2vcf)

    if os.path.isfile(join(vcfpath, (vcf[1] + '.somatic.enhanced.vcf'))):
        logger.info(f"VCF2VCF succeeded. {join(vcfpath, (vcf[1] + '.somatic.enhanced.vcf'))} exists!")
        subprocess.call(vcf2maf)

        if os.path.isfile(join(mafpath, (vcf[1] + '.somatic.maf'))):
            logger.info(f"VCF2MAF succeeded. {join(mafpath, (vcf[1] + '.somatic.maf'))} exists!")
            subprocess.call(maf2maf)

            if os.path.isfile(join(enhancedpath, (vcf[1] + '.somatic.enhanced.maf'))):
                logger.info(f"MAF2MAF succeeded. {join(mafpath, (vcf[1] + '.somatic.maf'))} exists!")
                success +=1

logger.info(f'Enhanced MAF files found: {success}')
for root, dirs, files in os.walk(enhancedpath):
    for file in files:
        logger.info(file)

print("done")

