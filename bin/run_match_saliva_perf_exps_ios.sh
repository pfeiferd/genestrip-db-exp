#!/bin/sh
set -e

scriptdir=$(dirname "$0")

cd $scriptdir/..

for t in ERR1395613 ERR1395610 SRR5571991 SRR5571990 SRR5571985;
  do
    rm -f ./data/projects/viral/csv/viral_match_${t}.csv
    /usr/bin/time -l mvn exec:exec@singlematch -Dname=viral -Dgoal=match -Dfqfile=${t}_1.fastq.gz  >& data/match_viral_${t}.log
  done
