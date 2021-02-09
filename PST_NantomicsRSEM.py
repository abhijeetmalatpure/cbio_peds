import os
from os.path import expanduser, join
import subprocess
from configparser import ConfigParser
import time, datetime
import logging
import socket
import pandas as pd
import concurrent.futures
import threading
import gzip
import glob
logger = None

# Create log directory
if not os.path.isdir('log'):
    os.mkdir('log')

# Before we do anything, set up logging.
epoch_time = time.mktime(datetime.datetime.now().timetuple())
log_filename = os.path.join(f"log/nantomics_rsem_{str(int(epoch_time))}.log")

print(f'Writing to {log_filename}')
logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s %(message)s", filename=log_filename)
logger = logging.getLogger()

print(f"Starting download on {socket.gethostname()} with process ID: {os.getpid()}")

# Read vcf metadata file into memory
vcf_metadata: pd.DataFrame = pd.read_csv('/N/project/phi_ri/cmg_dev_ri/germline_cmg_saliva_join_vcfs/mm_cmg_vcf_germline_files.csv', delimiter='|',
                           names=['vendor', 'sample_id', 'filename', 'patient_id', 'disease', 'uuid']).astype(str)

peds_metadata = pd.read_csv('/N/project/phi_ri/nantomics/peds_matching/Riley_Research_Report_UUIDs_sorted_by_contrast.tsv', delimiter='|', header = None)
peds_metadata.columns = peds_metadata.iloc[44]
peds_metadata['rsemFile'] = ''

for index, row in peds_metadata.iterrows():
    for name in glob.glob('/N/project/phi_ingest_nantomics/nantomics/nant/peds/pst-files-per-patient/*/BAM/*' + row['ReportUUID'] + '*rsem.txt.gz'):
        print(name)
        peds_metadata.iloc[index]['rsemFile'] = name


run = str(int(time.mktime(datetime.datetime.now().timetuple())))

rsem_root = '/N/project/phi_ingest_nantomics/nantomics/nant/peds/pst-files-per-patient'
sema = threading.Semaphore(value=5)
threads = list()

# Load samtools, bcftools
config = ConfigParser()
config.read('download.cfg')
libraries = dict(config.items('libraries'))

os.environ['LOADEDMODULES'] = libraries['loadedmodules']
os.environ['PATH'] = libraries['path']
os.environ['_LMFILES_'] = libraries['lmfiles']
os.environ['LD_LIBRARY_PATH'] = libraries['ld_library_path']

toolspath = join(expanduser("~"), "cbio_tools", "vcf2maf")
rsempath = "/N/project/phi_asha_archive/peds_pst/nantomics/tmp"
tmppath = "/N/project/phi_asha_archive/peds_pst/foundation/tmp"

[os.makedirs(path, exist_ok=True) for path in [rsempath]]

logger.info(f'Enhanced MAF files before processing: {len(os.listdir(rsempath))}')


def extract_rnaseq(vcf_filename, sample_id, vcf_type):
    print(vcf_filename)
    sema.acquire()
    tmpdir = join(tmppath, ("MG_" + sample_id.split('-')[0]), vcf_type)
    os.makedirs(tmpdir, exist_ok=True)
    with gzip.open(vcf_filename, mode='rt') as rsem:
        rnaseq = open(join(tmpdir, vcf_type + '.vcf'), 'wt')
        for line in rsem:
            rnaseq.write(line)
        rnaseq.close()


for index, row in vcf_metadata.iterrows():
    germ = threading.Thread(target=extract_rnaseq, args=(row.Germline, row.ContrastUUID, 'germline'))
    threads.append(germ)
    germ.start()

for thread in threads:
    threading.Thread.join(thread)


vcf_metadata.to_csv(os.path.join(rsempath, "vcf_metadata_" + run + ".csv"), sep="\t")

logger.info("PEDS Nantomics completed.")