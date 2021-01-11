import os
from os.path import expanduser, join, sep
import subprocess
from configparser import ConfigParser
import time, datetime
import logging
import socket
import pandas as pd
import gzip
import glob

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

print(f"Starting download on {socket.gethostname()} with process ID: {os.getpid()}")


config = ConfigParser()
config.read('download.cfg')
libraries = dict(config.items('libraries'))

os.environ['LOADEDMODULES'] = libraries['loadedmodules']
os.environ['PATH'] = libraries['path']
os.environ['_LMFILES_'] = libraries['lmfiles']
os.environ['LD_LIBRARY_PATH'] = libraries['ld_library_path']

toolspath = join(expanduser("~"), "cbio_tools", "vcf2maf")
vcfpath = "/N/slate/abhmalat/peds_pst/nantomics/vcf2mafConversion/vcf"
mafpath = "/N/slate/abhmalat/peds_pst/nantomics/vcf2mafConversion/maf"
enhancedpath = "/N/slate/abhmalat/peds_pst/nantomics/vcf2mafConversion/enhanced"

[os.makedirs(path, exist_ok=True) for path in [vcfpath, mafpath, enhancedpath]]


# Read the metadata file, find location of files using glob
vcf_metadata = pd.read_csv('/N/project/phi_ri/nantomics/peds_matching/Riley_Research_Report_UUIDs_sorted_by_contrast.tsv', delimiter='|', header = None)

vcf_metadata.columns = vcf_metadata.iloc[44]

vcf_metadata['Germline'] = ''
vcf_metadata['Somatic'] = ''

for index, row in vcf_metadata.iterrows():
    for name in glob.glob('/N/project/phi_ingest_nantomics/nantomics/nant/peds/pst-files-per-patient/*/VCF/*' + row['ContrastUUID'] + '*vcf.gz'):
        if name.endswith('germ.vcf.gz'):
            vcf_metadata.iloc[index]['Germline'] = name
        elif name.endswith('som.vcf.gz'):
            vcf_metadata.iloc[index]['Somatic'] = name

vcf_metadata.drop(44, axis=0, inplace=True)
vcf_metadata.reset_index(inplace=True, drop=True)


logger.info(f'Enhanced MAF files before processing: {len(os.listdir(enhancedpath))}')

run = str(int(time.mktime(datetime.datetime.now().timetuple())))
tmpdir = {'v2m': join(expanduser("~"), 'tmp', ("MG_" + '99b81dc8dd7d'), ("V2M_" + run))}
[os.makedirs(tmp, exist_ok=True) for tmp in tmpdir.values()]


def vcf2maf_call(vcf_filename, sample_id, vcf_type):
    with gzip.open(vcf_filename, mode='r') as vcf:
        vcf_extract = open(join(tmpdir['v2m'], 'germline.vcf'), 'wb')
        for line in vcf:
            vcf_extract.write(line)
        vcf_extract.close()

    id_type = '--normal-id' if vcf_type == 'germline' else '--tumor-id'

    maf_file = join(mafpath, sample_id + '.' + vcf_type + '.maf')
    vcf2maf = [
        'perl', join(toolspath, 'vcf2maf.pl'),
        '--input-vcf', vcf_extract.name,
        '--output-maf', maf_file,
        '--vep-path', join(toolspath, 'vep_hg19'),
        '--vep-data', join(expanduser("~"), '.vep'),
        '--vep-forks', '5',
        '--tmp-dir', tmpdir['v2m'],
        '--vcf-tumor-id', 'TUMOR',
        '--vcf-normal-id', 'NORMAL',
        id_type, sample_id,
        '--ncbi-build', 'GRCh37',
        '--ref-fasta',
        join(expanduser("~"),
             '.vep/homo_sapiens/102_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz')
    ]
    logger.info('VCF2MAF: ' + ' '.join(vcf2maf))

    os.unlink(vcf_extract.name)


# For each row unzip
for index, row in vcf_metadata:
    # Process germline file for sample
    vcf2maf_call(row.Germline, row.ContrastUUID, 'germline')
    vcf2maf_call(row.Somatic, row.ContrastUUID, 'somatic')


# evcffile = '/N/project/phi_ingest_nantomics/nantomics/nant/peds/pst-files-per-patient/' \
#            'c21b0ed9-2478-4ca1-a094-7dd29c773530_aQka3c2U3TGyMtL76aWGbJ/VCF/adaae6c5-5895-4ca3-a00f-34f68044524e_2017-09-04.germ.vcf.gz'
#
# id = evcffile.split('/')[-1].split('_')[0]
#
# maffile = join(mafpath, 'adaae6c5-5895-4ca3-a00f-34f68044524e_2017-09-04.germ.maf')
#
#
#     # file_content = vcf.read()
#     # fp.write(file_content)
#
#
#
# vcf2maf = [
#     'perl', join(toolspath, 'vcf2maf.pl'),
#     '--input-vcf', fp.name,
#     '--output-maf', maffile,
#     '--vep-path', join(toolspath, 'vep_hg19'),
#     '--vep-data', join(expanduser("~"), '.vep'),
#     '--vep-forks', '5',
#     '--tmp-dir', tmpdir['v2m'],
#     '--vcf-tumor-id', 'TUMOR',
#     '--vcf-normal-id', 'NORMAL',
#     '--normal-id', 'adaae6c5-5895-4ca3-a00f-34f68044524e',
#     '--ncbi-build', 'GRCh37',
#     '--ref-fasta',
#     join(expanduser("~"),
#          '.vep/homo_sapiens/102_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz')
# ]
#
#
# print(' '.join(vcf2maf))
