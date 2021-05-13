#!/usr/bin/env bash

# generate a picking list of files to download; can be passed as an argument ot transform.py

if [ $# -gt 0 ]; then
    OUTPUTSTEM="${1}"
    GUIDPATTERN="${1}*"
    GUIDPATTERNLEN=`expr length ${1}`

    echo "Selecting guids starting with $GUIDPATTERN (length of selection string = $GUIDPATTERNLEN )"

else
    echo "Selecting all guids"
    GUIDPATTERN="*"
    GUIDPATTERNLEN=0
    OUTPUTSTEM = "all"
fi

TARGETDIR="../testdata/"
SEARCHPATH="/mnt/microbio/ndm-hicf/ogre/pipeline_output/${GUIDPATTERN}/MAPPING/2e6b7bc7-f52c-4649-8538-c984ab3894bb_R00000039/STD/basecalls/${GUIDPATTERN}_v3.fasta.gz"
# get the file list
FILELIST="${OUTPUTSTEM}.txt"
rsync -r oxford.generic@158.119.199.95:${SEARCHPATH} > $FILELIST
NFILES=`cat $FILELIST | wc -l`
echo "The picking list contains $NFILES files in $FILELIST"

