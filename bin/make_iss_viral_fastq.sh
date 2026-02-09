#!/bin/sh
set -e

scriptdir=$(dirname "$0")

mkdir -p $scriptdir/..

cd $scriptdir/..

# Get viral data
mkdir -p iss
mkdir -p iss/viral
cp data/common/refseq/viral.1.1.genomic.fna.gz iss/viral
gunzip iss/viral/viral.1.1.genomic.fna.gz

iss generate --genomes iss/viral/viral.1.1.genomic.fna  --model miseq --output data/projects/viral/fastq/viral_iss_miseq_reads
iss generate --genomes iss/viral/viral.1.1.genomic.fna  --model hiseq --output data/projects/viral/fastq/viral_iss_hiseq_reads