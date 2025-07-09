#!/bin/sh
set -e

scriptdir=$(dirname "$0")

cd $scriptdir/..

file_path=./data/fastq

for t in ERR1395613 ERR1395610 SRR5571991 SRR5571990 SRR5571985;
  do
    ./cgmemtime/cgmemtime ./ku/krakenuniq/krakenuniq --report-file ${file_path}/ku_${t}.tsv --output off --threads 10 -db ./ku/viral_db ${file_path}/${t}_1.fastq.gz  >& data/match_ku_viral_${t}.log
  done
