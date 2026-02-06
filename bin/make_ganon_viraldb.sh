#!/bin/sh
set -e

scriptdir=$(dirname "$0")

mkdir -p $scriptdir/../ganon

cd $scriptdir/../ganon

# Create a viral database
mkdir -p viral_db

# Requires FTP to work.
ganon build --db-prefix viral_db --source refseq --organism-group viral --threads 24

# Run ganon on simulated virus fastq from paper:
ganon classify --db-prefix viral_db -s ../data/projects/viral/fastq/viral_fasta2fastq_fasta1.fastq.gz --output-all --output-all -o viralfastq1 --threads 32
