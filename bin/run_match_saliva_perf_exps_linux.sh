#!/bin/sh
set -e

scriptdir=$(dirname "$0")

cd $scriptdir/..

res_path=./results

for t in ERR1395613 ERR1395610 SRR5571991 SRR5571990 SRR5571985;
  do
    rm -f ./data/projects/viral/csv/viral_match_${t}.csv
    ./cgmemtime/cgmemtime mvn exec:exec@singlematch -Dname=viral -Dgoal=match -Dfqfile=${t}_1.fastq.gz  > ${res_path}/logs/match_viral_${t}.log
    mv ./data/projects/viral/csv/viral_match_${t}.csv ${res_path}/reports
  done
