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
from shutil import copyfile

# Create log directory
if not os.path.isdir('log'):
    os.mkdir('log')

# Before we do anything, set up logging.
epoch_time = time.mktime(datetime.datetime.now().timetuple())
log_filename = os.path.join(f"log/PST_nantomics_rsem_{str(int(epoch_time))}.log")

print(f'Writing to {log_filename}')
logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s %(message)s", filename=log_filename)
logger = logging.getLogger()

print(f"Process: {socket.gethostname()} with process ID: {os.getpid()}")

rsem_root = '/N/project/phi_ingest_nantomics/nantomics/nant/peds/pst-files-per-patient'
tmppath = "/N/project/phi_asha_archive/peds_pst/nantomics/tmp"
rsempath = "/N/project/phi_asha_archive/peds_pst/nantomics/rsem"
annotatedpath = "/N/project/phi_asha_archive/peds_pst/nantomics/rsem/annotated2"

[os.makedirs(path, exist_ok=True) for path in [tmppath, rsempath, annotatedpath]]

rsem_combined = "rsem_combined_2.csv"

run = str(int(time.mktime(datetime.datetime.now().timetuple())))

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

rsem_colnames = ["transcript","chromosome", "Strand", "transcript_start", "transcript_end",
                 "CDS_start", "CDS_end", "exon_number", "exon_start_positions", "exon_end_positions",
                 "gene_id", "transcript2", "protein_coding", "length", "effective_length", "expected_count",
                 "TPM", "FPKM", "IsoPct", "pme_expected_count", "pme_TPM", "pme_FPKM", "sample_id"]

def combine_rsem_files():
    combined_file = os.path.join(rsempath, rsem_combined)
    # if os.path.isfile(combined_file):
    #    return pd.read_csv(combined_file, delimiter=",")

    # Read vcf metadata file into memory
    vcf_metadata: pd.DataFrame = pd.read_csv('/N/project/phi_asha_archive/peds_pst/nantomics/vcf2mafConversion/maf/peds_vcf_metadata.csv', delimiter=',', header=0)

    peds_metadata = pd.read_csv('/N/project/phi_ri/nantomics/peds_matching/Riley_Research_Report_UUIDs_sorted_by_contrast.tsv', delimiter='|', header = None)
    peds_metadata.columns = peds_metadata.iloc[44]
    vcf_metadata['rsem'] = ''

    for index, row in vcf_metadata.iterrows():
        if isinstance(row.somatic_tumor_rna, str):
            for name in glob.glob('/N/project/phi_ingest_nantomics/nantomics/nant/peds/pst-files-per-patient/*/BAM/*' + row['somatic_tumor_rna'] + '*.rsem.txt.gz'):
                # print(name)
                vcf_metadata.loc[index, 'rsem'] = name

    vcf_metadata.dropna(axis=0, subset=['somatic_tumor_rna'], inplace=True)

    [os.makedirs(path, exist_ok=True) for path in [rsempath]]

    # vcf_metadata['rsem_filename'] = vcf_metadata['rsem'].str.split('/')[-1]

    nant_rsem = pd.DataFrame(data={}, columns=rsem_colnames)

    for index, row in vcf_metadata.iterrows():
        fname = row.rsem.split('/')[-1][0:-3]
        row.rsem_filename = os.path.join(rsempath, fname)
        temp_rsem = pd.read_csv(row.rsem_filename, sep="\t")
        temp_rsem.columns = nant_rsem.columns[:-1]

        #print(temp_rsem.head(10))

        # Save a copy with
        temp_rsem.to_csv(os.path.join(annotatedpath, fname),
                         sep="\t", header=True, index=False)

        temp_rsem["sample_id"] = row.tumor_id
        nant_rsem = nant_rsem.append(temp_rsem, ignore_index=True)

        if not os.path.isfile(os.path.join(rsempath, fname)):
            copyfile(row.rsem, os.path.join(rsempath, fname))


    return nant_rsem

nant_combined_rsem = combine_rsem_files()

#nant_combined_rsem.to_csv(os.path.join(rsempath, "combined.rsem"), index=False, header=True, sep="\t")
logger.info("PEDS Nantomics RSEM completed.")