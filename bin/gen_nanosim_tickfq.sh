#!/bin/sh
set -e

scriptdir=$(dirname "$0")

cd $scriptdir/..
basedir=$(pwd)

mkdir -p nanosim/viral
mkdir -p nanosim/tick-borne

# Does not work -c nanosim/viral/training, -c does not work at all
read_analysis.py metagenome --fastq -gl data/projects/viral/csv/viral_nanosim.tsv -i data/fastq/tick1.fastq.gz  -t 24
read_analysis.py metagenome --fastq -gl data/projects/tick-borne/csv/tick-borne_nanosim.tsv -i data/fastq/tick1.fastq.gz  -t 24

# simulator.py metagenome --fastq -gl data/projects/viral/csv/viral_nanosim.tsv -t 24 -a abundance.tsv
