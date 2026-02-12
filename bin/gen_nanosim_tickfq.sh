#!/bin/sh
set -e

scriptdir=$(dirname "$0")

cd $scriptdir/..
basedir=$(pwd)

nanosimdir = ${basedir}/../Nanosim/src

mkdir -p nanosim/viral
mkdir -p nanosim/tick-borne

conda activate nanosim

# Does not work -c nanosim/viral/training, -c does not work at all
${nanosimdir}/read_analysis.py metagenome -q --fastq -gl data/projects/viral/csv/viral_nanosim.tsv -i data/fastq/tick1.fastq.gz  -t 24
#${nanosimdir}/read_analysis.py quantify -e meta -gl data/projects/viral/csv/viral_nanosim.tsv -i data/fastq/tick1.fastq.gz  -t 24
#${nanosimdir}/read_analysis.py metagenome --fastq -gl data/projects/tick-borne/csv/tick-borne_nanosim.tsv -i data/fastq/tick1.fastq.gz  -t 24

reads=$(zcat tick1.fastq.gz | wc -l | awk '{print $1/4}')
# TODO generate abundanc.tsv from training_quantification.tsv using ${reads}
${nanosimdir}/simulator.py metagenome --fastq -gl data/projects/viral/csv/viral_nanosim.tsv -t 24 -a abundance.tsv
mv simulated_sample0_unaligned_reads.fastq nanosim/viral/sim_tick1.fastq