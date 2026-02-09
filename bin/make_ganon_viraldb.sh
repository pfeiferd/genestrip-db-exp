#!/bin/sh
set -e

scriptdir=$(dirname "$0")

mkdir -p $scriptdir/..

cd $scriptdir/..

# Create a viral database
mkdir -p ../ganon/viral_db
mkdir -p ../ganon/viral_lowfp_db

# Requires FTP to work.
#ganon build --db-prefix viral_db --source refseq --organism-group viral --threads 24

# Build ganon database
ganon build-custom --input-file ../data/projects/viral/csv/viral_ganon.tsv --taxonomy-files ../data/common/nodes.dmp ../data/common/names.dmp --db-prefix ../ganon/viral_db --level species --threads 32
ganon build-custom --input-file ../data/projects/viral/csv/viral_ganon.tsv --taxonomy-files ../data/common/nodes.dmp ../data/common/names.dmp --db-prefix ../ganon/viral_lowfp_db --level species --threads 32 --max-fp 0.0000001

# Run ganon on simulated virus fastq from paper:
ganon classify --db-prefix ../ganon/viral_db -s ../data/projects/viral/fastq/viral_fasta2fastq_fasta1.fastq.gz --output-all --output-all -o ../ganon/viral_fastq1 --threads 32
ganon classify --db-prefix ../ganon/viral_db -s ../data/projects/viral/fastq/viral_iss_hiseq_reads_R1.fastq ../data/projects/viral/fastq/viral_iss_hiseq_reads_R2.fastq --output-all --output-all -o ../ganon/viral_iss_hiseq --threads 32
ganon classify --db-prefix ../ganon/viral_db -s ../data/projects/viral/fastq/viral_iss_miseq_viral_reads_R1.fastq ../data/projects/viral/fastq/viral_iss_miseq_reads_R2.fastq --output-all --output-all -o ../ganon/viral_iss_miseq --threads 32

ganon classify --db-prefix ../ganon/viral_lowfp_db -s ../data/projects/viral/fastq/viral_fasta2fastq_fasta1.fastq.gz --output-all --output-all -o ../ganon/viral_lowfp_fastq1 --threads 32
ganon classify --db-prefix ../ganon/viral_lowfp_db -s ../data/projects/viral/fastq/viral_iss_hiseq_reads_R1.fastq ../data/projects/viral/fastq/viral_iss_hiseq_reads_R2.fastq --output-all --output-all -o ../ganon/viral_lowfp_iss_hiseq --threads 32
ganon classify --db-prefix ../ganon/viral_lowfp_db -s ../data/projects/viral/fastq/viral_iss_miseq_viral_reads_R1.fastq ../data/projects/viral/fastq/viral_iss_miseq_reads_R2.fastq --output-all --output-all -o ../ganon/viral_lowfp_iss_miseq --threads 32
