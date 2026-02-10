#!/bin/sh
set -e

scriptdir=$(dirname "$0")

cd $scriptdir/..
basedir=$(pwd)

## Ganon ##

# Run ganon on simulated virus fastq from paper:
for db in viral viral_lowfp human_virus human_virus_lowfp;
  do
    ganon classify --db-prefix ./ganon/${db}_db -s ./data/fastq/viral_fasta2fastq_fasta1.fastq.gz --output-all -o ./ganon/${db}_fastq1 --threads 32
    ganon classify --db-prefix ./ganon/${db}_db -s ./data/fastq/viral_iss_hiseq_reads_R1.fastq ./data/fastq/viral_iss_hiseq_reads_R2.fastq --output-all --output-all -o ./ganon/${db}_iss_hiseq --threads 32
    ganon classify --db-prefix ./ganon/${db}_db -s ./data/fastq/viral_iss_miseq_reads_R1.fastq ./data/fastq/viral_iss_miseq_reads_R2.fastq --output-all --output-all -o ./ganon/${db}_iss_miseq --threads 32

    for id in ERR1395613 ERR1395610 SRR5571991 SRR5571990 SRR5571985;
    do
      ganon classify --db-prefix ./ganon/${db}_db -s ./data/fastq/${id}_1.fastq.gz ./data/fastq/${id}_2.fastq.gz --output-all --output-all -o ./ganon/${db}_${id} --threads 32
    done
  done

## KrakenUniq ##

for db in viral human_virus;
  do
    ./ku/krakenuniq/krakenuniq --threads 10 -db ./ku/${db}_db ./data/fastq/viral_fasta2fastq_fasta1.fastq.gz --output ./ku/${db}_fastq1.tsv
    ./ku/krakenuniq/krakenuniq --threads 10 -db ./ku/${db}_db ./data/fastq/viral_iss_hiseq_reads_R1.fastq ./data/fastq/viral_iss_hiseq_reads_R2.fastq --output ./ku/${db}_iss_hiseq.tsv
    ./ku/krakenuniq/krakenuniq --threads 10 -db ./ku/${db}_db ./data/fastq/viral_iss_miseq_reads_R1.fastq ./data/fastq/viral_iss_miseq_reads_R2.fastq --output ./ku/${db}_iss_miseq.tsv

    for id in ERR1395613 ERR1395610 SRR5571991 SRR5571990 SRR5571985;
    do
      ./ku/krakenuniq/krakenuniq --only-classified-output --threads 10 -db ./ku/${db}_db ./data/fastq/${id}_1.fastq ./data/fastq/${id}_2.fastq --output ./ku/${db}_${id}.tsv
    done
  done

