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
import threading

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
vcfpath = "/N/project/phi_asha_archive/peds_pst/nantomics/vcf2mafConversion/vcf"
mafpath = "/N/project/phi_asha_archive/peds_pst/nantomics/vcf2mafConversion/maf"
enhancedpath = "/N/project/phi_asha_archive/peds_pst/nantomics/vcf2mafConversion/enhanced"
tmppath = "/N/project/phi_asha_archive/peds_pst/nantomics/vcf2mafConversion/tmp"

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
vcf_metadata['normal_id'] = ""
vcf_metadata['tumor_id'] = ""

logger.info(f'MAF files generated: {len(os.listdir(mafpath))}')

run = str(int(time.mktime(datetime.datetime.now().timetuple())))

sema = threading.Semaphore(value=5)
threads = list()

def vcf2maf_call(vcf_filename, sample_id, vcf_type):
    print(vcf_filename)
    sema.acquire()
    tmpdir = join(tmppath, ("MG_" + sample_id.split('-')[0]), vcf_type)
    os.makedirs(tmpdir, exist_ok=True)

    with gzip.open(vcf_filename, mode='rt') as vcf:
        vcf_extract = open(join(tmpdir, vcf_type + '.vcf'), 'wt')
        for line in vcf:
            vcf_extract.write(line)
            if line.startswith('##SAMPLE=<ID=NORMAL,'):
                # print(line)
                normal_id = dict(item.split('=') for item in line[line.find('<') + 1:line.find('>')].replace("\"", '').split(','))['File'].split('.')[0]
            if line.startswith('##SAMPLE=<ID=TUMOR,'):
                # print(line)
                tumor_id = dict(item.split('=') for item in line[line.find('<') + 1:line.find('>')].replace("\"", '').split(','))['File'].split('.')[0]

        if vcf_type == 'somatic':
            vcf_metadata.loc[vcf_metadata.Somatic == vcf_filename, 'tumor_id'] = tumor_id
        else:
            vcf_metadata.loc[vcf_metadata.Germline == vcf_filename, 'normal_id'] = normal_id
        vcf_extract.close()

    # maf_file = join(mafpath, sample_id + '.' + vcf_type + '.maf')
    # vcf2maf = [
    #     'perl', join(toolspath, 'vcf2maf.pl'),
    #     '--input-vcf', vcf_extract.name,
    #     '--output-maf', maf_file,
    #     '--vep-path', join(toolspath, 'vep_hg19'),
    #     '--vep-data', join(expanduser("~"), '.vep'),
    #     '--vep-forks', '5',
    #     '--tmp-dir', tmpdir,
    #     '--vcf-tumor-id', 'TUMOR',
    #     '--vcf-normal-id', 'NORMAL',
    #     '--normal-id', normal_id,
    #     '--tumor-id', tumor_id,
    #     '--ncbi-build', 'GRCh37',
    #     '--ref-fasta',
    #     join(expanduser("~"),
    #          '.vep/homo_sapiens/102_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz')
    # ]

    # if not os.path.isfile(maf_file):
    #     logger.info(f"{maf_file} did not exist. Running vcf2maf to generate it.")
    #     logger.info('VCF2MAF: ' + ' '.join(vcf2maf))
    #     subprocess.call(vcf2maf, stdout=open(os.devnull, 'w'))
    #
    #     if(os.path.isfile(maf_file)):
    #         logger.info(f"Success for {vcf_extract.name}!")
    #     else:
    #         logger.info(f"Failure for {vcf_extract.name}!")
    # else:
    #     logger.info(f"{maf_file} exists. Moving on.")

    # Delete the uncompressed file
    if os.path.isfile(vcf_extract.name):
        for item in os.listdir(tmpdir):
            if item.endswith(".vcf"):
                logger.info(f'Deleting {os.path.join(tmpdir, item)}')
                if os.path.isfile(os.path.join(tmpdir, item)):
                    os.unlink(os.path.join(tmpdir, item))
                else:
                    logger.info(f'File {os.path.join(tmpdir, item)} does not exist')

    sema.release()


logger.info(f"Total Samples: {vcf_metadata.shape[0]}")

# For each row unzip, run vcf2maf, delete temp file
for index, row in vcf_metadata.iterrows():
    germ = threading.Thread(target=vcf2maf_call, args=(row.Germline, row.ContrastUUID, 'germline'))
    som = threading.Thread(target=vcf2maf_call, args=(row.Somatic, row.ContrastUUID, 'somatic'))
    threads.append(germ)
    threads.append(som)
    germ.start()
    som.start()

for thread in threads:
    threading.Thread.join(thread)

vcf_metadata.to_csv(os.path.join(mafpath, "vcf_metadata_" + run + ".csv"), sep="\t")

logger.info("PEDS Nantomics completed.")
