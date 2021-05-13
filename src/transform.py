""" generates mixture marked fasta files from a 'picklist' made by make_picklist.sh """
import os
import io
import gzip
import sys
import pathlib
import logging
import collections
from vcfScan import FastaMixtureMarker, vcfScan
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna

try:
    filelist=sys.argv[1]
except IndexError:
    print("Needs one command line argument, the picking list to use")
    exit(1)
if os.path.exists(filelist):
    print("Using picking list:",filelist)
else:
    print("No valid picking list provided: file does not exist")
    exit(1)

TARGETDIR="/srv/data/mixfiles/mfasta/"

# make sure a target directory exists
pathlib.Path(TARGETDIR).mkdir(parents=True, exist_ok=True)

# run a check to see how many samples need processing.
missing= set()
with open(filelist,'r') as f:
    for line in f.readlines():
        GUID =line[46:(46+36)] 
        fasta_output_filename = os.path.join(TARGETDIR, '{0}.mfasta.gz'.format(GUID))
        if not os.path.exists(fasta_output_filename):
              missing.add(GUID)
print("There are {0} files found which are not processed".format(len(missing)))


print('setting up vcfScan object..')
# use the baseCounts4 tag to identify high quality bases
# don't report minor variants present at less that 5% frequency (done here simply to speed up computations)
v = vcfScan(expectedErrorRate = 0.001, infotag = 'BaseCounts4', report_minimum_maf = 0.05, compute_pvalue = False)     

# we define one region for each base of the genome
for i in range(4411532):     
    v.add_roi(str(1+i),set([1+i]))

nSkipped = 0
nFailed = 0
nSucceeded = 0 
failures = []
for GUID in missing:
        FASTAFILE=os.path.join(TARGETDIR, "{0}.fasta.gz".format(GUID))
        VCFFILE=os.path.join(TARGETDIR, "{0}.vcf.gz".format(GUID))
        fasta_output_filename = os.path.join(TARGETDIR, '{0}.mfasta.gz'.format(GUID))
        cmd1="scp oxford.generic@158.119.199.95:/mnt/microbio/ndm-hicf/ogre/pipeline_output/{0}/MAPPING/2e6b7bc7-f52c-4649-8538-c984ab3894bb_R00000039/STD/basecalls/{0}_v3.fasta.gz {1}".format(GUID,FASTAFILE)
        cmd2="scp oxford.generic@158.119.199.95:/mnt/microbio/ndm-hicf/ogre/pipeline_output/{0}/MAPPING/2e6b7bc7-f52c-4649-8538-c984ab3894bb_R00000039/STD/basecalls/{0}_v3.vcf.gz {1}".format(GUID,VCFFILE)

        fasta_success =  os.path.exists(FASTAFILE)
        if not fasta_success:
            os.system(cmd1)
        composition_good = None
        fasta_success =  os.path.exists(FASTAFILE)
        if fasta_success:
            # run a composition check on the fasta file.
            try:
                
                with gzip.open(FASTAFILE, "rt") as handle:
                    for input_record in SeqIO.parse(handle, "fasta"):
                        composition = collections.Counter(list(input_record.seq))
                        # check whether all four bases are called
                        expected_bases = set(['A','C','G','T'])
                        observed_bases = set(composition.keys()).intersection(expected_bases)
                        composition_good =  (len(observed_bases)==4)
                        print(GUID, composition, composition_good)

            except OSError:
                # may no be a gzipped file
                print("Read failed with OSError for {0}".format(FASTAFILE))
                os.unlink(FASTAFILE)
                     
            if composition_good is True:		# analyse and annotate the fasta filename
                    print("Composition check is good")
                    vcf_success = os.path.exists(VCFFILE)
                    if not vcf_success: 
                       os.system(cmd2)

                    vcf_success = os.path.exists(VCFFILE)

                    if os.path.exists(fasta_output_filename):
                        nSkipped+=1
                        print("Output file exists")
                    else:
                        print(GUID,"being processed")
                        if not os.path.exists(FASTAFILE):
                           os.system(cmd1)
                        if not os.path.exists(VCFFILE):
                           os.system(cmd2)

                        if os.path.exists(VCFFILE) and os.path.exists(FASTAFILE):
                            print('parsing vcf')
                            if v.parse(vcffile = VCFFILE):		# succeeded
                                mixfile = os.path.join(TARGETDIR,'{0}.txt'.format(GUID))
                                v.bases.to_csv(mixfile, index=None)

                                # the fasta file used contains high confidence single base base calls.  Example is from PHE TB pipeline https://github.com/oxfordmmm/CompassCompact
                                fmm = FastaMixtureMarker(expectedErrorRate = 0.001, mlp_cutoff=6.65, clustering_cutoff = 10, min_maf=0)
                                df, seq = fmm.mark_mixed(FASTAFILE, mixfile)

                                if df is not None:      # it succeeded
                                    iupac = ['r','R','w','W','y','Y','m','M','s','S','k','K']
                                    nMixed=0
                                    for item in iupac:
                                        nMixed = nMixed + seq.count(item)
                                    print("There were {0} mixed bases".format(nMixed))

                                    # write fasta
                                    record = SeqRecord(Seq(seq,
                                               generic_dna),
                                               id=GUID, name="{0}|Consensus with mixtures marked".format(GUID),
                                               description="consensus sequence with mixed bases recorded as iupac codes")
                                    # write to string   
                                    with io.StringIO("") as f:
                                        SeqIO.write(record, f, 'fasta')
                                        fout = f.getvalue()
                                        with gzip.open(fasta_output_filename, mode='wb') as fg:
                                           fg.write(fout.encode('utf-8'))
                                    print("Wrote mfasta to {0}".format(fasta_output_filename))

                                    print("Consensus fasta file produced")
                                    nSucceeded +=1
                                else:
                                    print("Fastamixturemasker returned None.  setting composition to false")
                                    composition_good= False
                            else:
                                composition_good = False
                                print("VCF parse failed; composition marked as bad")
                                
                                
                          
            if composition_good is False:       # composition is not good; sequence has no data, or the vcf file cannot be parsed
                    # write fasta
                    print("Composition check failed, or VCF parse failed")
                    record = SeqRecord(input_record.seq,
                               id=GUID, name="{0}|Consensus".format(GUID),
                               description="consensus sequence is pipeline derived")
                    # write to string   
                    with io.StringIO("") as f:
                        SeqIO.write(record, f, 'fasta')
                        fout = f.getvalue()
                        with gzip.open(fasta_output_filename, mode='wb') as fg:
                           fg.write(fout.encode('utf-8'))
                    print("Used pipeline consensus.  Wrote fasta -> mfasta to {0}".format(fasta_output_filename))
            if os.path.exists(VCFFILE):
                print("Removing vcf file")
                os.unlink(VCFFILE)

            if os.path.exists(FASTAFILE):
                print("Keeping original fasta file")
                #os.unlink(FASTAFILE)
    
print("Complete.  Skipped {0} files, which had already been processed.".format(nSkipped))
print("Successfully processed {0} files.".format(nSucceeded))
print("Failed to process {0} files, which had already been processed.".format(nFailed))

