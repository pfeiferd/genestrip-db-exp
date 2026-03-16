#!/bin/sh
set -e

scriptdir=$(dirname "$0")

cd $scriptdir/..
basedir=$(pwd)

# Use 'nanosimdir' to point to the python scripts for NanoSim.
nanosimdir=${nanosim:-${basedir}/../NanoSim/src}

#conda init
#conda activate nanosim

if [ -e ./data/projects/tick-borne/fasta/GCF_016920785.2_ASM1692078v2_genomic.fna ]
then
    echo Fna file exists # Do nothing
else
  cp ./data/common/fasta/GCF_016920785.2_ASM1692078v2_genomic.fna.gz ./data/projects/tick-borne/fasta
  gunzip ./data/projects/tick-borne/fasta/GCF_016920785.2_ASM1692078v2_genomic.fna.gz
fi

# This does not work: Abundances become 99% tick DNA - we don't want that.
#sed -i '1i\6945\t./data/projects/tick-borne/fasta/GCF_016920785.2_ASM1692078v2_genomic.fna' data/projects/tick-borne/csv/tick-borne_nanosim.tsv

# Generate fastq files with tick DNA only
for file in tick1 tick2 tick3 tick4 tick5 tick6 tick7 tick8
do
  ${nanosimdir}/read_analysis.py metagenome -q --fastq -gl data/projects/tick-borne/csv/tick_only_nanosim.tsv -i data/fastq/${file}.fastq.gz  -t 24

  reads=1000000 #$(zcat data/fastq/${file}.fastq.gz | wc -l | awk '{print $1/4}')
  sed -i "s/Abundance/$reads/g" training_quantification.tsv
  ${nanosimdir}/simulator.py metagenome --seed 42 --fastq -gl data/projects/tick-borne/csv/tick_only_nanosim.tsv -t 24 -a training_quantification.tsv
  mv simulated_sample0_aligned_reads.fastq data/fastq/${file}_tick_only_sim.fastq
  rm training*
  rm reference_metagenome.fasta
done

# Mixing in tick DNA - too complicated (for now) - we leave it as it is.
for x in ; do
# Generate fastq files with bacterial DNA only
for file in tick1 tick2 tick3 tick4 tick5 tick6 tick7 tick8
do
  ${nanosimdir}/read_analysis.py metagenome -q --fastq -gl data/projects/tick-borne/csv/tick-borne_nanosim.tsv -i data/fastq/${file}.fastq.gz  -t 24

  reads=1000000 #$(zcat data/fastq/${file}.fastq.gz | wc -l | awk '{print $1/4}')
  sed -i "s/Abundance/$reads/g" training_quantification.tsv
  ${nanosimdir}/simulator.py metagenome --seed 42 --fastq -gl data/projects/tick-borne/csv/tick-borne_nanosim.tsv -t 24 -a training_quantification.tsv
  mv simulated_sample0_aligned_reads.fastq data/fastq/${file}_sim.fastq
  rm training*
  rm reference_metagenome.fasta
done

# Generate fastq with mixed bacterial DNA and tick DNA
# (We could not get NanoSim to do this any other way because tick DNA got 99% abundance always.)
for file in tick1 tick2 tick3 tick4 tick5 tick6 tick7 tick8
do
  cat data/fastq/${file}_sim.fastq data/fastq/${file}_tick_only_sim.fastq > data/fastq/${file}_mixed_sim.fastq
done

done