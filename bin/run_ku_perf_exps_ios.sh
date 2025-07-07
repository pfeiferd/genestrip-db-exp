#!/bin/sh
set -e

scriptdir=$(dirname "$0")

cd $scriptdir/..

ps=24G
file_path=./data/fastq

for t in tick1 tick2 tick3 tick4 tick5 tick6 tick7 tick8;
  do
    /usr/bin/time -l ./ku/krakenuniq/krakenuniq --threads 9 --preload-size ${ps} -db ./ku/microbial_db ${file_path}/${tb}.fastq.gz  >& data/match_ku_mb_${t}.log
  done

cd ${file_path}

rm -f ${file_path}/allticks.fastq.gz
cat tick1.fastq.gz tick2.fastq.gz tick3.fastq.gz tick4.fastq.gz tick5.fastq.gz tick6.fastq.gz tick7.fastq.gz tick8.fastq.gz > allticks.fastq.gz

/usr/bin/time -l ./ku/krakenuniq/krakenuniq --threads 9 --preload-size ${ps} -db ./ku/microbial_db ${file_path}/allticks.fastq.gz >& data/match_ku_mb_allticks.log
