#!/bin/sh
set -e

scriptdir=$(dirname "$0")

cd $scriptdir/..
basedir=$(pwd)

nanosimdir=${basedir}/../NanoSim/src

conda init
conda activate nanosim

for file in tick1 tick2
do
  ${nanosimdir}/read_analysis.py metagenome -q --fastq -gl data/projects/viral/csv/viral_nanosim.tsv -i data/fastq/${file}.fastq.gz  -t 24
  #${nanosimdir}/read_analysis.py metagenome --fastq -gl data/projects/tick-borne/csv/tick-borne_nanosim.tsv -i data/fastq/tick1.fastq.gz  -t 24

  reads=$(zcat data/fastq/${file}.fastq.gz | wc -l | awk '{print $1/4}')
  sed -i "s/Abundance/$reads/g" training_quantification.tsv
  ${nanosimdir}/simulator.py metagenome --fastq -gl data/projects/viral/csv/viral_nanosim.tsv -t 24 -a training_quantification.tsv
  mv simulated_sample0_unaligned_reads.fastq data/fastq/${file}_sim.fastq
done