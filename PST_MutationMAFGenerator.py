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

#Before we do anything, set up logging.
epoch_time = time.mktime(datetime.datetime.now().timetuple())
log_filename = os.path.join(f"log/mafconverter_{str(int(epoch_time))}.log")

print(f'Writing to {log_filename}')
logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s %(message)s", filename=log_filename)
logger = logging.getLogger()

print(f"Starting download on {socket.gethostname()} with process ID: {os.getpid()}")
# Each list within vcf_files has-
# 1. File path
# 2. Sample name to be used - "filename"
# 3. Sample name found in vcf file
# The output maf and enhanced maf will have names
# filename.somatic.enhanced.vcf, filename.somatic.maf and filename.somatic.enhanced.maf (or germline)
vcfs = []
header = []
# Uncomment this and remove the other os.walk for processing all files
# for root, dirs, files in os.walk("/N/slate/abhmalat/peds_pst/ashion/files/subjects"):
#         for file in files:
#             if file.endswith(".vcf") and root.split(sep)[-1] == 'normalized':
#                 vcftype = 'somatic' if file.endswith('somatic.vcf') else 'germline'
#                 with open(join(root, file), 'r') as content:
#                     for line in content:
#                         line = line.strip()
#                         if line.startswith("#CHROM\t"):
#                             vcfs.append([join(root, file),
#                                          file.split('.')[0],
#                                          line.split("\t")[-1],
#                                          vcftype])
#                             break

foi = ['C051_0040_032766_T1_K1ID2_ps20201022180611.germline.vcf',
'C051_0037_031894_T1_K1ID2_ps20201014112210.germline.vcf',
'C051_0043_035165_T1_K1ID2_ps20201129101455.germline.vcf',
'C051_0008_020553_T1_K1ID2_ps20200318083449.germline.vcf',
'C051_0028_027595_T1_K1ID2_ps20200813161941.germline.vcf',
'C051_0006_020025_T1_K1ID2_ps20200303145923.germline.vcf',
'C051_0021_025563_T1_K1ID2_ps20200707145222.germline.vcf',
'C051_0016_022944_T1_K1ID2_ps20200521163846.germline.vcf',
'C051_0032_029353_T1_K1ID2_ps20200909130840.germline.vcf',
'C051_0033_030037_T1_K1ID2_ps20200915123511.germline.vcf',
'C051_0012_021539_T1_K1ID2_ps20200407144223.germline.vcf',
'C051_0017_022784_T1_K1ID2_ps20200519142241.germline.vcf',
'C051_0022_025621_T1_K1ID2_ps20200710144706.germline.vcf',
'C051_0027_027303_T1_K1ID2_ps20200805153358.germline.vcf',
'C051_0005_020208_T1_K1ID2_ps20200310140712.germline.vcf',
'C051_0019_023207_T1_K1ID2_ps20200528134022.germline.vcf',
'C051_0030_029040_T1_K1ID2_ps20200903134757.germline.vcf',
'C051_0004_018940_T1_K1ID2_ps20200211195123.germline.vcf',
'C051_0029_028371_T1_K1ID2_ps20200821153024.germline.vcf',
'C051_0036_031751_T1_K1ID2_ps20201009131438.germline.vcf',
'C051_0015_022721_T1_K1ID2_ps20200519142239.germline.vcf',
'C051_0003_018942_T1_K1ID2_ps20200211195116.germline.vcf',
'C051_0010_020979_T1_K1ID2_ps20200325090600.germline.vcf',
'C051_0041_034401_T1_K1ID2_ps20201114064729.germline.vcf',
'C051_0013_021790_T1_K1ID2_ps20200421094355.germline.vcf']

for root, dirs, files in os.walk("/N/slate/abhmalat/peds_pst/ashion/files/subjects"):
        for file in files:
            if file in foi and root.split(sep)[-1] == 'normalized':
                vcftype = 'somatic' if file.endswith('somatic.vcf') else 'germline'
                with open(join(root, file), 'r') as content:
                    for line in content:
                        line = line.strip()
                        if line.startswith("#CHROM\t"):
                            vcfs.append([join(root, file),
                                         file.split('.')[0],
                                         line.split("\t")[-1],
                                         vcftype])
                            break


#logger.info(f'Normalized VCF files found:{len(vcfs)} ')
logger.info(f'Normalized VCF files found:{len(vcfs)} ')
#logger.info(vcfs)

[os.makedirs(join(expanduser("~"), 'tmp', (vcf[1] + vcf[3])),
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


logger.info(f'Enhanced MAF files before processing: {len(os.listdir(enhancedpath))}')

success = 0
for vcf in vcfs:
    # Remove temp files generated in the previous run
    tmpdir = join(expanduser("~"), 'tmp', (vcf[1] + vcf[3]))
    for tmpfile in os.listdir(tmpdir):
        if os.path.isfile(tmpfile) or os.path.islink(tmpfile):
            os.unlink(tmpfile)
            logger.info(f"File {tmpfile} deleted from {tmpdir}")

    evcffile = join(vcfpath, (vcf[1] + '.' + vcf[3] + '.enhanced.vcf'))
    maffile = join(mafpath, (vcf[1] + '.' + vcf[3] + '.maf'))
    enhancedfile = join(enhancedpath, (vcf[1] + '.' + vcf[3] + '.enhanced.maf'))

    if os.path.isfile(evcffile):
        os.unlink(evcffile)
    if os.path.isfile(maffile):
        os.unlink(maffile)
    if os.path.isfile(enhancedfile):
        os.unlink(enhancedfile)

    vcf2vcf = [
        'perl', join(toolspath, 'vcf2vcf.pl'),
        '--input-vcf', vcf[0],
        '--output-vcf', evcffile,
        '--vcf-tumor-id' if vcf[3] == 'somatic' else '--vcf-normal-id', vcf[2],
        '--new-tumor-id' if vcf[3] == 'somatic' else '--new-normal-id', vcf[1],
        '--ref-fasta',
        join(expanduser("~"),
             '.vep/homo_sapiens/102_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz')]

    vcf2maf = [
            'perl', join(toolspath, 'vcf2maf.pl'),
            '--input-vcf', evcffile,
            '--output-maf', maffile,
            '--vep-path', join(toolspath, 'vep_hg19'),
            '--vep-data', join(expanduser("~"), '.vep'),
            '--vep-forks', '5',
            '--tmp-dir', tmpdir,
            '--vcf-tumor-id' if vcf[3] == 'somatic' else '--vcf-normal-id', vcf[1],
            '--tumor-id' if vcf[3] == 'somatic' else '--normal-id', vcf[1],
            '--ncbi-build', 'GRCh37',
            '--ref-fasta',
            join(expanduser("~"),
                 '.vep/homo_sapiens/102_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz')
        ]

    maf2maf = [
                'perl', join(toolspath, 'maf2maf.pl'),
                '--input-maf', maffile,
                '--output-maf', enhancedfile,
                '--vep-path', join(toolspath, 'vep_hg19'),
                '--vep-data', join(expanduser("~"), '.vep'),
                '--vep-forks', '5',
                '--tmp-dir', tmpdir,
                '--ncbi-build', 'GRCh37',
                '--ref-fasta',
                join(expanduser("~"),
                     '.vep/homo_sapiens/102_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz')
            ]

    logger.info("VCF2VCF: " + ' '.join(vcf2vcf))
    logger.info("VCF2MAF: " + ' '.join(vcf2maf))
    logger.info("MAF2MAF: " + ' '.join(maf2maf))

    subprocess.call(vcf2vcf)

    if os.path.isfile(evcffile):
        logger.info(f"VCF2VCF succeeded. {evcffile} exists!")
        subprocess.call(vcf2maf)

        if os.path.isfile(maffile):
            logger.info(f"VCF2MAF succeeded. {maffile} exists!")
            subprocess.call(maf2maf)

            if os.path.isfile(enhancedfile):
                logger.info(f"MAF2MAF succeeded. {enhancedfile} exists!")
                success +=1
            else:
                logger.info(
                    f"MAF2MAF FAILED FOR {maffile}")
        else:
            logger.info(
                f"VCF2MAF FAILED FOR {evcffile}")
    else:
        logger.info(
            f"VCF2VCF FAILED FOR {vcf[0]}")


logger.info(f'Enhanced MAF files found: {len(os.listdir(enhancedpath))}. Success: {success}')# {success}')
for root, dirs, files in os.walk(enhancedpath):
    for file in files:
        logger.info(file)

print("done")

