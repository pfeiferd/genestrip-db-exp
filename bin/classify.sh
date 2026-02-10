#!/bin/sh
set -e

scriptdir=$(dirname "$0")

cd $scriptdir/..
basedir=$(pwd)

## Ganon ##

# Run ganon on simulated virus fastq from paper:
ganon classify --db-prefix ../ganon/viral_db -s ../data/projects/viral/fastq/viral_fasta2fastq_fasta1.fastq.gz --output-all --output-all -o ../ganon/viral_fastq1 --threads 32
ganon classify --db-prefix ../ganon/viral_db -s ../data/projects/viral/fastq/viral_iss_hiseq_reads_R1.fastq ../data/projects/viral/fastq/viral_iss_hiseq_reads_R2.fastq --output-all --output-all -o ../ganon/viral_iss_hiseq --threads 32
ganon classify --db-prefix ../ganon/viral_db -s ../data/projects/viral/fastq/viral_iss_miseq_reads_R1.fastq ../data/projects/viral/fastq/viral_iss_miseq_reads_R2.fastq --output-all --output-all -o ../ganon/viral_iss_miseq --threads 32

ganon classify --db-prefix ../ganon/viral_lowfp_db -s ../data/projects/viral/fastq/viral_fasta2fastq_fasta1.fastq.gz --output-all --output-all -o ../ganon/viral_lowfp_fastq1 --threads 32
ganon classify --db-prefix ../ganon/viral_lowfp_db -s ../data/projects/viral/fastq/viral_iss_hiseq_reads_R1.fastq ../data/projects/viral/fastq/viral_iss_hiseq_reads_R2.fastq --output-all --output-all -o ../ganon/viral_lowfp_iss_hiseq --threads 32
ganon classify --db-prefix ../ganon/viral_lowfp_db -s ../data/projects/viral/fastq/viral_iss_miseq_reads_R1.fastq ../data/projects/viral/fastq/viral_iss_miseq_reads_R2.fastq --output-all --output-all -o ../ganon/viral_lowfp_iss_miseq --threads 32

ganon classify --db-prefix ../ganon/human_virus_db -s ../data/projects/viral/fastq/viral_fasta2fastq_fasta1.fastq.gz --output-all --output-all -o ../ganon/human_virus_fastq1 --threads 32
ganon classify --db-prefix ../ganon/human_virus_db -s ../data/projects/viral/fastq/viral_iss_hiseq_reads_R1.fastq ../data/projects/viral/fastq/viral_iss_hiseq_reads_R2.fastq --output-all --output-all -o ../ganon/human_virus_iss_hiseq --threads 32
ganon classify --db-prefix ../ganon/human_virus_db -s ../data/projects/viral/fastq/viral_iss_miseq_reads_R1.fastq ../data/projects/viral/fastq/viral_iss_miseq_reads_R2.fastq --output-all --output-all -o ../ganon/human_virus_iss_miseq --threads 32

ganon classify --db-prefix ../ganon/human_virus_lowfp_db -s ../data/projects/viral/fastq/viral_fasta2fastq_fasta1.fastq.gz --output-all --output-all -o ../ganon/human_virus_lowfp_fastq1 --threads 32
ganon classify --db-prefix ../ganon/human_virus_lowfp_db -s ../data/projects/viral/fastq/viral_iss_hiseq_reads_R1.fastq ../data/projects/viral/fastq/viral_iss_hiseq_reads_R2.fastq --output-all --output-all -o ../ganon/human_virus_lowfp_iss_hiseq --threads 32
ganon classify --db-prefix ../ganon/human_virus_lowfp_db -s ../data/projects/viral/fastq/viral_iss_miseq_reads_R1.fastq ../data/projects/viral/fastq/viral_iss_miseq_reads_R2.fastq --output-all --output-all -o ../ganon/human_virus_lowfp_iss_miseq --threads 32
