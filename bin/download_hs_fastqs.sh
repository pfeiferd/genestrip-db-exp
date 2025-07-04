#!/bin/sh
set -e

scriptdir=$(dirname "$0")

cd $scriptdir/..

bin=./sra-tools/sratoolkit/bin
file_path=./data/fastq

# ERX1462737
# ERX1462740
# SRX2830683
# SRX2830684
# SRX2830689
for id in SRR5571991 SRR5571990 SRR5571985; # ERR1395610 ERR1395613
  do
      $bin/prefetch $id --max-size 200g -O $file_path
      $bin/fasterq-dump $file_path/$id -O $file_path
      gzip $file_path/${id}_?.fastq
      rm -rf $file_path/${id}
  done
