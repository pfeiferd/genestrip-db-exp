#!/bin/sh
set -e

scriptdir=$(dirname "$0")

mkdir -p $scriptdir/..

cd $scriptdir/..

# Create a viral database
mkdir -p ../ganon/viral_db

# Requires FTP to work.
#ganon build --db-prefix viral_db --source refseq --organism-group viral --threads 24

# Build ganon database
ganon build-custom --input-file ../data/projects/viral/csv/viral_ganon.tsv --taxonomy-files ../data/common/nodes.dmp ../data/common/names.dmp --db-prefix ../ganon/viral_db --level species --threads 32

# Run ganon on simulated virus fastq from paper:
ganon classify --db-prefix ../ganon/viral_db -s ../data/projects/viral/fastq/viral_fasta2fastq_fasta1.fastq.gz --output-all --output-all -o ../ganon/viralfastq1 --threads 32
ganon classify --db-prefix ../ganon/viral_db -s ../data/projects/viral/fastq/iss_hiseq_viral_reads_R1.fastq ../data/projects/viral/fastq/iss_hiseq_viral_reads_R2.fastq --output-all --output-all -o ../ganon/iss_hiseq_viral --threads 32
ganon classify --db-prefix ../ganon/viral_db -s ../data/projects/viral/fastq/iss_miseq_viral_reads_R1.fastq ../data/projects/viral/fastq/iss_miseq_viral_reads_R2.fastq --output-all --output-all -o ../ganon/iss_miseq_viral --threads 32
