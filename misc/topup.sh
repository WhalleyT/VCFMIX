#!/usr/bin/env bash

# don't run if this script is already running - disabled
# we are using flock to effect this, within the calling script.
#dupl=$(ps -x |  grep "bash ./topup.sh" | grep -v grep | wc -l)
#echo -e "There are $dupl related processes already running"
#if [ "$dupl" -gt 0 ]; then
#  echo -e "script running already; terminating";
#  exit 0
#fi

# truncate log files to 1M ## TODO
#tail /srv/data/mixfiles/log/*.log

# first, get a list of all guids processed by the pipeline.
# this is  relatively quick process
echo `date;who`
echo `date;echo "topup.sh | finding target files"`

# get guids beginning with 0,1,2, etc
cd /home/phe.gov.uk/david.wyllie/VCFMIX/src/
rm *.txt		# remove all output files

./picklist.sh 0
./picklist.sh 1
./picklist.sh 2
./picklist.sh 3
./picklist.sh 4
./picklist.sh 5
./picklist.sh 6
./picklist.sh 7
./picklist.sh 8
./picklist.sh 9
./picklist.sh a
./picklist.sh b
./picklist.sh c
./picklist.sh d
./picklist.sh e
./picklist.sh f
cat *.txt > all.txt 		# concatenate into one big file

# ensure these files are converted to .mfasta.gz (mixtures marked) format
echo `date;echo "topup.sh | transforming target files"`

python3 /home/phe.gov.uk/david.wyllie/VCFMIX/src/transform.py all.txt

# push any outstanding items into findNeighbour4, keeping the server up to date
echo `date;echo "topup.sh | Push into server"`

cd /home/phe.gov.uk/david.wyllie/findNeighbour4/src

# when run with cron, need to specify full path.
/home/phe.gov.uk/david.wyllie/.local/bin/pipenv run python3 demo_pheall.py

echo `date;echo "topup.sh | Finished"`

