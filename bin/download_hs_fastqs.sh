#!/bin/sh
set -e

scriptdir=$(dirname "$0")

cd $scriptdir/..

bin=./sra-tools/sratoolkit/bin
file_path=../../data/fastq

# ERX1462737
$bin/prefetch ERR1395610 --max-size 200g -O $file_path
# ERX1462740
$bin/prefetch ERR1395613 --max-size 200g -O $file_path
# SRX2830683
$bin/prefetch SRR5571991 --max-size 200g -O $file_path
# SRX2830684
$bin/prefetch SRR5571990 --max-size 200g -O $file_path
# SRX2830689
$bin/prefetch SRR5571985 --max-size 200g -O $file_path

$bin/fasterq-dump $file_path/ERR1395610 -O $file_path
$bin/fasterq-dump $file_path/ERR1395613 -O $file_path
$bin/fasterq-dump $file_path/SRR5571991 -O $file_path
$bin/fasterq-dump $file_path/SRR5571990 -O $file_path
$bin/fasterq-dump $file_path/SRR5571985 -O $file_path
