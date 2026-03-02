#!/bin/sh
set -e

scriptdir=$(dirname "$0")

cd $scriptdir/..
basedir=$(pwd)

nanosimdir=${basedir}/../NanoSim/src

#conda init
#conda activate nanosim

# cp ./data/common/fasta/GCF_016920785.2_ASM1692078v2_genomic.fna.gz ./data/projects/tick-borne/fasta
# gunzip ./data/projects/tick-borne/fasta/GCF_016920785.2_ASM1692078v2_genomic.fna.gz
# sed -i '1i\6945\t./data/projects/tick-borne/fasta/GCF_016920785.2_ASM1692078v2_genomic.fna' data/projects/tick-borne/csv/tick-borne_nanosim.tsv

for file in tick1 tick2 tick3 tick4 tick5 tick6 tick7 tick8
do
  #${nanosimdir}/read_analysis.py metagenome -q --fastq -gl data/projects/viral/csv/viral_nanosim.tsv -i data/fastq/${file}.fastq.gz  -t 24
  ${nanosimdir}/read_analysis.py metagenome -q --fastq -gl data/projects/tick-borne/csv/tick-borne_nanosim.tsv -i data/fastq/${file}.fastq.gz  -t 24

  reads=$(zcat data/fastq/${file}.fastq.gz | wc -l | awk '{print $1/4}')
  sed -i "s/Abundance/$reads/g" training_quantification.tsv
  #${nanosimdir}/simulator.py metagenome --fastq -gl data/projects/viral/csv/viral_nanosim.tsv -t 24 -a training_quantification.tsv
  ${nanosimdir}/simulator.py metagenome --fastq -gl data/projects/tick-borne/csv/tick-borne_nanosim.tsv -t 24 -a training_quantification.tsv
  mv simulated_sample0_aligned_reads.fastq data/fastq/${file}_sim.fastq
  rm training*
  rm reference_metagenome.fasta
done