#!/bin/sh
set -e

scriptdir=$(dirname "$0")

cd $scriptdir/..

# Generate simulated viral data over all RefSeq viruses using iss
mkdir -p iss
mkdir -p iss/viral
cp data/common/refseq/viral.1.1.genomic.fna.gz iss/viral
gunzip iss/viral/viral.1.1.genomic.fna.gz

iss generate --genomes iss/viral/viral.1.1.genomic.fna  --model miseq --output data/fastq/viral_iss_miseq_reads
iss generate --genomes iss/viral/viral.1.1.genomic.fna  --model hiseq --output data/fastq/viral_iss_hiseq_reads

# Download ONT-based tick related fastqs using wget
file_path=./data/fastq

cd $file_path

wget https://www.be-md.ncbi.nlm.nih.gov/Traces/sra-reads-be/fastq?acc=SRR17281117 -O tick1.fastq.gz
wget https://www.be-md.ncbi.nlm.nih.gov/Traces/sra-reads-be/fastq?acc=SRR17281105 -O tick2.fastq.gz
wget https://www.be-md.ncbi.nlm.nih.gov/Traces/sra-reads-be/fastq?acc=SRR17281103 -O tick3.fastq.gz
wget https://www.be-md.ncbi.nlm.nih.gov/Traces/sra-reads-be/fastq?acc=SRR17281101 -O tick4.fastq.gz
wget https://www.be-md.ncbi.nlm.nih.gov/Traces/sra-reads-be/fastq?acc=SRR17281100 -O tick5.fastq.gz
wget https://www.be-md.ncbi.nlm.nih.gov/Traces/sra-reads-be/fastq?acc=SRR17281099 -O tick6.fastq.gz
wget https://www.be-md.ncbi.nlm.nih.gov/Traces/sra-reads-be/fastq?acc=SRR17281116 -O tick7.fastq.gz
wget https://www.be-md.ncbi.nlm.nih.gov/Traces/sra-reads-be/fastq?acc=SRR17281115 -O tick8.fastq.gz


# Download human saliva fastq using sra-tools
bin=./sra-tools/sratoolkit/bin
file_path=./data/fastq

# ERX1462737
# ERX1462740
# SRX2830683
# SRX2830684
# SRX2830689
for id in ERR1395613 ERR1395610 SRR5571991 SRR5571990 SRR5571985;
  do
      $bin/prefetch $id --max-size 200g -O $file_path
      $bin/fasterq-dump $file_path/$id -O $file_path
      gzip $file_path/${id}_?.fastq
      rm -rf $file_path/${id}
  done

