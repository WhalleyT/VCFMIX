""" Example of use of FastaMixtureMarker """
import os
import sys
import urllib.request
from pathlib import Path
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna
SOURCE_DIR = Path(__file__).parent.absolute()
src_dir = os.path.join(SOURCE_DIR, '..', 'src')
testdata_dir = os.path.join(SOURCE_DIR, '..', 'data', 'testdata')
sys.path.append(src_dir)
from vcfScan import FastaMixtureMarker, vcfScan

print('set up vcfScan object.. (only needs to be done once)')
# use the baseCounts4 tag to identify high quality bases
# don't report minor variants present at less that 5% frequency (done here simply to speed up computations)
v = vcfScan(expectedErrorRate=0.001, infotag='BaseCounts4', report_minimum_maf=0.05, compute_pvalue=False)

# we define one region for each base of the genome
for i in range(4411532):
    v.add_roi(str(1 + i), set([1 + i]))

# retrieve an input file; parse input file
test_vcf_file = os.path.join(testdata_dir, '52858be2-7020-4b7f-acb4-95e00019a7d7_v3.vcf.gz')
if not os.path.exists(test_vcf_file):
    url = 'https://ora.ox.ac.uk/objects/uuid:5e4ec1f8-e212-47db-8910-161a303a0757/download_file?file_format=x-tar&safe_filename=52858be2-7020-4b7f-acb4-95e00019a7d7_v3.vcf.gz&type_of_work=Dataset'
    urllib.request.urlretrieve(url, test_vcf_file)

res = v.parse(vcffile=test_vcf_file)
print("Parse complete; writing output")

# make sure a target directory exists
targetdir = os.path.join(SOURCE_DIR, 'unitTest_tmp')  # a writeable directory
Path(targetdir).mkdir(parents=True, exist_ok=True)

# write mixed bases to a csv file
mixfile = os.path.join(targetdir, 'output_table.txt')
v.bases.to_csv(mixfile, index=None)

# the fasta file used contains high confidence single base base calls. Example is from PHE TB pipeline https://github.com/oxfordmmm/CompassCompact
fastafile = os.path.join(testdata_dir, '52858be2-7020-4b7f-acb4-95e00019a7d7_v3.fasta')
fmm = FastaMixtureMarker(expectedErrorRate=0.001, mlp_cutoff=6.65, clustering_cutoff=10, min_maf=0)
df, seq = fmm.mark_mixed(fastafile, mixfile)

iupac = ['A', 'C', 'G', 'T', 'r', 'R', 'w', 'W', 'y', 'Y', 'm', 'M', 's', 'S', 'k', 'K']
resDict = {}
for item in iupac:
    resDict[item] = seq.count(item)
    print("There were {1} mixed bases of type {0}".format(item, resDict[item]))

# write fasta
record = SeqRecord(Seq(seq, generic_dna), id="test_id", name="fmm_output_test",
                    description="consensus sequence with mixed bases recorded as iupac codes")

fasta_output_filename = os.path.join(targetdir, 'test_output.fasta')

with open(fasta_output_filename, 'w') as f:
    SeqIO.write(record, f, 'fasta')
