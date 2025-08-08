#!/bin/bash
set -e

scriptdir=$(dirname "$0")

cd $scriptdir/..

file_path=./data/fastq
res_path=./results

for t in ERR1395613 ERR1395610 SRR5571991 SRR5571990 SRR5571985;
  do
    ./cgmemtime/cgmemtime ./ku/krakenuniq/krakenuniq --report-file ${res_path}/reports/ku_${t}.tsv --output off --threads 10 -db ./ku/viral_db ${file_path}/${t}_1.fastq.gz  >& ${res_path}/logs/match_ku_viral_${t}.log
  done
