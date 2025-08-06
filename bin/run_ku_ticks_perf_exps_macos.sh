#!/bin/sh
set -e

scriptdir=$(dirname "$0")

cd $scriptdir/..

res_path=./results

ps=24G
file_path=./data/fastq
# Beware: 10 Thread is correct here as Genestrip is configured with 9 consumer thread + 1 reading thread...
for t in tick1 tick2 tick3 tick4 tick5 tick6 tick7 tick8;
  do
    /usr/bin/time -l ./ku/krakenuniq/krakenuniq --report-file ${res_path}/reports/ku_${t}.tsv --output off --threads 10 --preload-size ${ps} -db ./ku/microbial_db ${file_path}/${t}.fastq.gz  >& ${res_path}/logs/match_ku_mb_${t}.log
  done

/usr/bin/time -l ./ku/krakenuniq/krakenuniq --report-file ${res_path}/reports/ku_allticks.tsv --output off --threads 10 --preload-size ${ps} -db ./ku/microbial_db ${file_path}/tick1.fastq.gz ${file_path}/tick2.fastq.gz ${file_path}/tick3.fastq.gz ${file_path}/tick4.fastq.gz ${file_path}/tick5.fastq.gz ${file_path}/tick6.fastq.gz ${file_path}/tick7.fastq.gz ${file_path}/tick8.fastq.gz >& ${res_path}/logs/match_ku_mb_allticks.log
