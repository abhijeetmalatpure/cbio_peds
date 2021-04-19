import os
from os.path import expanduser, join, sep
import subprocess
from configparser import ConfigParser
import time, datetime
import logging
import socket
import pandas as pd

logger = None

# Create log directory
if not os.path.isdir('log'):
    os.mkdir('log')

# Before we do anything, set up logging.
epoch_time = time.mktime(datetime.datetime.now().timetuple())
log_filename = os.path.join(f"log/PSTFoundationMafConverter_{str(int(epoch_time))}.log")

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
vcf_metadata = pd.DataFrame(columns=['filename', 'sample_id', 'vcf_tumor', 'vcftype'])
vcfs = []
header = []
# Uncomment this and remove the other os.walk for processing all files
for root, dirs, files in os.walk("/N/slate/abhmalat/peds_pst/foundation/files/foundation"):
    for file in files:
        if file.endswith(".vcf"):
            vcftype = 'somatic'
            with open(join(root, file), 'r') as content:
                for line in content:
                    line = line.strip()
                    if line.startswith("#CHROM\t"):
                        filename = join(root, file)
                        tumor_id = file.split('.')[0]
                        file_sample = line.split("\t")[-1]
                        vcf_metadata = vcf_metadata.append(pd.Series([filename, tumor_id, file_sample, vcftype],
                                                                     index=vcf_metadata.columns), ignore_index=True)
                        break

logger.info(f'VCF files found:{len(vcf_metadata)}')



# Load samtools, bcftools
config = ConfigParser()
config.read('download.cfg')
libraries = dict(config.items('libraries'))

os.environ['LOADEDMODULES'] = libraries['loadedmodules']
os.environ['PATH'] = libraries['path']
os.environ['_LMFILES_'] = libraries['lmfiles']
os.environ['LD_LIBRARY_PATH'] = libraries['ld_library_path']

toolspath = join(expanduser("~"), "cbio_tools", "vcf2maf")
mafpath = "/N/project/phi_asha_archive/peds_pst/foundation/maf2/somatic"
#enhancedpath = "/N/slate/abhmalat/peds_pst/foundation/vcf2mafConversion/enhanced"
tmppath = "/N/project/phi_asha_archive/peds_pst/foundation/tmp"

[os.makedirs(path, exist_ok=True) for path in [mafpath]]

logger.info(f'MAF files before processing: {len(os.listdir(mafpath))}')

run = str(int(time.mktime(datetime.datetime.now().timetuple())))
vcf_metadata.to_csv(os.path.join(mafpath, "vcf_metadata_foundation_" + run + ".csv"), sep="\t")


success = 0
# for vcf in vcfs:
for index, vcf in vcf_metadata.iterrows():
    logger.info(f"{index}. Starting with {vcf.filename}")
    run = str(int(time.mktime(datetime.datetime.now().timetuple())))
    tmpdir = join(tmppath, ("MG_" + vcf.sample_id))

    # Create tmp directories
    os.makedirs(tmpdir, exist_ok=True)

    maffile = join(mafpath, (vcf.sample_id + '.' + vcf.vcftype + '.maf'))

    if os.path.isfile(maffile):
        logger.info(f"{vcf.filename} already exists. Skipping")
        continue

    vcf2maf = [
        'perl', join(toolspath, 'vcf2maf.pl'),
        '--input-vcf', vcf.filename,
        '--output-maf', maffile,
        '--vep-path', join(toolspath, 'vep'),
        '--vep-data', join(expanduser("~"), '.vep'),
        '--vep-forks', '5',
        '--tmp-dir', tmpdir,
        '--vcf-tumor-id', vcf.vcf_tumor,
        '--tumor-id', vcf.sample_id,
        '--custom-enst', '/N/u/abhmalat/Carbonate/cbio_tools/vcf2maf/data/isoform_overrides_uniprot',
        '--remap-chain', join(toolspath, 'data', 'hg19_to_hg38.chain'),
        '--ref-fasta',
        join(expanduser("~"),
             '.vep/homo_sapiens/broad_hg38/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta')
    ]

    logger.info("VCF2MAF: " + ' '.join(vcf2maf))

    subprocess.call(vcf2maf)

    if os.path.isfile(maffile):
        logger.info(f"VCF2MAF succeeded. {maffile} exists!")
        success += 1
    else:
        logger.info(f"VCF2MAF FAILED FOR {vcf.filename}")

    for f in os.listdir(tmpdir):
        if os.path.isfile(f) and f.endswith('.vcf'):
            logger.info(f"Removing {f}")
            os.unlink(f)


logger.info(f'MAF files found: {len(os.listdir(mafpath))}. Success: {success}')  # {success}')
for root, dirs, files in os.walk(join(mafpath, 'somatic')):
    for file in files:
        logger.info(file)

print("Done with Foundation somatic files")
