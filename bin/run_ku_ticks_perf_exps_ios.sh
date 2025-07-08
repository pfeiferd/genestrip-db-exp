#!/bin/sh
set -e

scriptdir=$(dirname "$0")

cd $scriptdir/..

ps=16G
file_path=./data/fastq

for t in tick1 tick2 tick3 tick4 tick5 tick6 tick7 tick8;
  do
    /usr/bin/time -l ./ku/krakenuniq/krakenuniq --report-file ${file_path}/ku_${t}.tsv --output off --threads 9 --preload-size ${ps} -db ./ku/microbial_db ${file_path}/${t}.fastq.gz  >& data/match_ku_mb_${t}.log
  done

cd ${file_path}

/usr/bin/time -l ./ku/krakenuniq/krakenuniq --report-file ${file_path}/ku_allticks.tsv --output off --threads 9 --preload-size ${ps} -db ./ku/microbial_db ${file_path}/tick1.fastq.gz ${file_path}/tick2.fastq.gz ${file_path}/tick3.fastq.gz ${file_path}/tick4.fastq.gz ${file_path}/tick5.fastq.gz ${file_path}/tick6.fastq.gz ${file_path}/tick7.fastq.gz ${file_path}/tick8.fastq.gz >& data/match_ku_mb_allticks.log
