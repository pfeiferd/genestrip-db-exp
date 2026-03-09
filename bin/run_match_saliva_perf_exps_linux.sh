#!/bin/sh
set -e

scriptdir=$(dirname "$0")

cd $scriptdir/..
basedir=$(pwd)

res_path=${basedir}/results

# Genestrip
for t in ERR1395613 ERR1395610 SRR5571991 SRR5571990 SRR5571985;
  do
    rm -f ./data/projects/viral/csv/viral_match_${t}.csv
    ./cgmemtime/cgmemtime mvn exec:exec@singlematch -Dname=viral -Dgoal=match -Dfqfile=${t}_1.fastq.gz  > ${res_path}/logs/match_viral_${t}.log
    mv ./data/projects/viral/csv/viral_match_${t}.csv ${res_path}/reports
  done

# KrakenUniq
file_path=./data/fastq

for t in ERR1395613 ERR1395610 SRR5571991 SRR5571990 SRR5571985;
  do
    # Suppress a lot of output to be fair
    ./cgmemtime/cgmemtime ./ku/krakenuniq/krakenuniq --report-file ./ku/viral_perf_ku_${t}.tsv --output off --threads 10 -db ./ku/viral_db ${file_path}/${t}_1.fastq.gz  > ${res_path}/logs/match_ku_viral_${t}.log
  done

# Kraken 2
for t in ERR1395613 ERR1395610 SRR5571991 SRR5571990 SRR5571985;
  do
    # Suppress a lot of output to be fair
    ./cgmemtime/cgmemtime ./k2/kraken2/k2 classify --threads 10 --db ./k2/viral_db --output - --report ./k2/viral_perf_k2_${t}.tsv ${file_path}/${t}_1.fastq.gz  > ${res_path}/logs/match_k2_viral_${t}.log
  done

# Ganon
for t in ERR1395613 ERR1395610 SRR5571991 SRR5571990 SRR5571985;
  do
    # Suppress a lot of output to be fair
    ./cgmemtime/cgmemtime ganon classify --db-prefix ./ganon/viral_db -o ./ganon/viral_perf_ganon_${t} --threads 10 -s ./data/fastq/${t}_1.fastq.gz > ${res_path}/logs/match_ganon_viral_${t}.log
  done

