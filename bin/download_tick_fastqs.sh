#!/bin/sh
set -e

scriptdir=$(dirname "$0")

cd $scriptdir/..

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
