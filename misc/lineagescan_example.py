import os
from src.vcfScan import lineageScan

# create a lineagescan object;
v = lineageScan()

# scan an input file
inputfile=os.path.join('..', 'data', 'testdata','52858be2-7020-4b7f-acb4-95e00019a7d7_v3.vcf.gz')
res = v.parse(vcffile = inputfile, guid='528')

# print details of the regions scanned
print(v.regions_stats)

# export details of the regions scanned
outputfile = os.path.join('unitTest_tmp','528.txt')
v.region_stats.to_csv(outputfile)

# compute F2 and F47 statistics (see publication)
summary1 = v.f_statistics()