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
vcf_metadata['Somatic'] = ""

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

evcffile = '/N/project/phi_ingest_nantomics/nantomics/nant/peds/pst-files-per-patient/' \
           '0159ccb0-bde8-446f-b241-79643ca5185d_hXqJfQ9jHH7LXu2BA87wv3/VCF/8792c312-e3cc-4e96-b0ec-99b81dc8dd7d_2018-06-30.germ.vcf.gz'

maffile = join(mafpath, ('8792c312-e3cc-4e96-b0ec-99b81dc8dd7d_2018-06-30.germ.vcf.gz' + '.maf'))

with gzip.open(evcffile, mode='r') as vcf:
    file_content = vcf.read()
    fp = open(join(tmpdir['v2m'], 'tmp.vcf'), 'wb')
    fp.write(file_content)



vcf2maf = [
    'perl', join(toolspath, 'vcf2maf.pl'),
    '--input-vcf', fp.name,
    '--output-maf', maffile,
    '--vep-path', join(toolspath, 'vep_hg19'),
    '--vep-data', join(expanduser("~"), '.vep'),
    '--vep-forks', '5',
    '--tmp-dir', tmpdir['v2m'],
    '--vcf-tumor-id', 'TUMOR',
    '--vcf-normal-id', 'NORMAL',
    '--tumor-id', 'TUMOR',
    '--normal-id', '99b81dc8dd7d',
    '--ncbi-build', 'GRCh37',
    '--ref-fasta',
    join(expanduser("~"),
         '.vep/homo_sapiens/102_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz')
]

print("VCF2MAF: " + ' '.join(vcf2maf))
