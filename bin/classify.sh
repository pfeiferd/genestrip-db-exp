#!/bin/sh
set -e

scriptdir=$(dirname "$0")

cd $scriptdir/..
basedir=$(pwd)

## Genestrip ##

#mvn exec:exec@matchrep -Dname=viral -Dgoal=match -Dfqmap=viral_acc_comp.txt
#mvn exec:exec@matchrep -Dname=human_virus -Dgoal=match -Dfqmap=viral_acc_comp.txt
#mvn exec:exec@matchrep2 -Dname=viral -Dgoal=match -Dfqmap=saliva.txt
#mvn exec:exec@matchrep2 -Dname=human_virus -Dgoal=match -Dfqmap=saliva.txt

## Ganon ##

# Run ganon on simulated virus fastq from paper:
for db in viral viral_lowfp human_virus human_virus_lowfp;
  do
    ganon classify --db-prefix ./ganon/${db}_db -s ./data/fastq/viral_fasta2fastq_fasta1.fastq.gz --output-all -o ./ganon/${db}_fastq1 --threads 32
    ganon classify --db-prefix ./ganon/${db}_db -s ./data/fastq/viral_iss_hiseq_reads_R1.fastq ./data/fastq/viral_iss_hiseq_reads_R2.fastq --output-all --output-all -o ./ganon/${db}_iss_hiseq --threads 32
    ganon classify --db-prefix ./ganon/${db}_db -s ./data/fastq/viral_iss_miseq_reads_R1.fastq ./data/fastq/viral_iss_miseq_reads_R2.fastq --output-all --output-all -o ./ganon/${db}_iss_miseq --threads 32
  done

# Now low fp for human saliva - we don't need it and it is so slow.
for db in viral human_virus;
  do
    for id in ERR1395613 # ERR1395610 SRR5571991 SRR5571990 SRR5571985;
    do
      ganon classify --db-prefix ./ganon/${db}_db -s ./data/fastq/${id}_1.fastq.gz ./data/fastq/${id}_2.fastq.gz --output-all --output-all -o ./ganon/${db}_${id} --threads 32
    done
  done

## KrakenUniq ##

for db in viral human_virus;
  do
    ./ku/krakenuniq/krakenuniq --threads 10 --db ./ku/${db}_db --output ./ku/${db}_fastq1.tsv ./data/fastq/viral_fasta2fastq_fasta1.fastq.gz
    ./ku/krakenuniq/krakenuniq --threads 10 --db ./ku/${db}_db --output ./ku/${db}_iss_hiseq.tsv ./data/fastq/viral_iss_hiseq_reads_R1.fastq ./data/fastq/viral_iss_hiseq_reads_R2.fastq
    ./ku/krakenuniq/krakenuniq --threads 10 --db ./ku/${db}_db --output ./ku/${db}_iss_miseq.tsv ./data/fastq/viral_iss_miseq_reads_R1.fastq ./data/fastq/viral_iss_miseq_reads_R2.fastq

   for id in ERR1395613 # ERR1395610 SRR5571991 SRR5571990 SRR5571985;
    do
      ./ku/krakenuniq/krakenuniq --only-classified-output --threads 10 --db ./ku/${db}_db --output ./ku/${db}_${id}.tsv ./data/fastq/${id}_1.fastq.gz ./data/fastq/${id}_2.fastq.gz
    done
  done

## Kraken2 ##

for db in viral human_virus;
  do
    ./k2/kraken2/k2 classify --threads 10 --db ./k2/${db}_db --output ./k2/${db}_fastq1.tsv ./data/fastq/viral_fasta2fastq_fasta1.fastq.gz
    ./k2/kraken2/k2 classify --threads 10 --db ./k2/${db}_db --output ./k2/${db}_iss_hiseq.tsv ./data/fastq/viral_iss_hiseq_reads_R1.fastq ./data/fastq/viral_iss_hiseq_reads_R2.fastq
    ./k2/kraken2/k2 classify --threads 10 --db ./k2/${db}_db --output ./k2/${db}_iss_miseq.tsv ./data/fastq/viral_iss_miseq_reads_R1.fastq ./data/fastq/viral_iss_miseq_reads_R2.fastq

   for id in ERR1395613 # ERR1395610 SRR5571991 SRR5571990 SRR5571985;
    do
      # First restrict to classified fastq file, then generate output to reduce effective output size.
      ./k2/kraken2/k2 classify --threads 10 --db ./k2/${db}_db --classified-out ./k2/classified_${id}.fastq --output - ./data/fastq/${id}_1.fastq.gz ./data/fastq/${id}_2.fastq.gz
      ./k2/kraken2/k2 classify --threads 10 --db ./k2/${db}_db ./k2/classified_${id}.fastq --output ./k2/${db}_${id}.tsv
    done
  done
