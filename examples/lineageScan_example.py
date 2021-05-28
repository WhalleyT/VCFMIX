""" Example of use of lineageScan """
import os
import sys
import urllib.request
import json
from pathlib import Path
SOURCE_DIR = Path(__file__).parent.absolute()
vcfmix_dir = os.path.join(SOURCE_DIR, '..', 'vcfmix')
testdata_dir = os.path.join(SOURCE_DIR, '..', 'data', 'testdata')
sys.path.append(vcfmix_dir)
from vcfScan import lineageScan

# create a lineagescan object;
v = lineageScan()

# retrieve an input file; parse input file
test_vcf_file = os.path.join(testdata_dir, '52858be2-7020-4b7f-acb4-95e00019a7d7_v3.vcf.gz')
if not os.path.exists(test_vcf_file):
    url = 'https://ora.ox.ac.uk/objects/uuid:5e4ec1f8-e212-47db-8910-161a303a0757/download_file?file_format=x-tar&safe_filename=52858be2-7020-4b7f-acb4-95e00019a7d7_v3.vcf.gz&type_of_work=Dataset'
    urllib.request.urlretrieve(url, test_vcf_file)

res = v.parse(vcffile=test_vcf_file, sample_id='528')

# print details of the regions scanned
print(v.region_stats)

# export details of the regions scanned
outputfile = os.path.join(SOURCE_DIR, 'examples_output', '528.txt')
v.region_stats.to_csv(outputfile)

# compute F2 and F47 statistics (see publication)
summary1 = v.f_statistics()
print(summary1)

# save F2 and F47 statistics to json
f_stats_outputfile = os.path.join(SOURCE_DIR, 'examples_output', '528_f_stats.json')
with open(f_stats_outputfile, 'w') as f: 
    json.dump(summary1, f)
