#!/bin/sh
set -e

scriptdir=$(dirname "$0")

cd $scriptdir/..
basedir=$(pwd)

res_path=${basedir}/results

for x in ;
do

# Genestrip
for t in tick1 tick2 tick3 tick4 tick5 tick6 tick7 tick8;
  do
    rm -f ./data/projects/tick-borne/csv/tick-borne_match_${t}.csv
    ./cgmemtime/cgmemtime mvn exec:exec@singlematch -Dname=tick-borne -Dgoal=match -Dfqfile=${t}.fastq.gz  > ${res_path}/logs/match_g_tick-borne_${t}.log
    mv ./data/projects/tick-borne/csv/tick-borne_match_${t}.csv ${res_path}/reports
  done

./cgmemtime/cgmemtime mvn exec:exec@match -Dname=tick-borne -Dgoal=match -Dfqmap=ticks.txt  > ${res_path}/logs/match_g_tick-borne_allticks.log

# KrakenUniq
file_path=${basedir}/data/fastq
ps=96G
# Beware: 10 Thread is correct here as Genestrip is configured with 9 consumer thread + 1 reading thread...
for t in tick1 tick2 tick3 tick4 tick5 tick6 tick7 tick8;
  do
    ./cgmemtime/cgmemtime ./ku/krakenuniq/krakenuniq --report-file ./ku/tick_perf_ku_${t}.tsv --output off --threads 10 --preload-size ${ps} -db ./ku/microbial_db ${file_path}/${t}.fastq.gz  > ${res_path}/logs/match_ku_mb_${t}.log
  done

./cgmemtime/cgmemtime ./ku/krakenuniq/krakenuniq --report-file ./ku/tick_perf_ku_allticks.tsv --output off --threads 10 --preload-size ${ps} -db ./ku/microbial_db ${file_path}/tick1.fastq.gz ${file_path}/tick2.fastq.gz ${file_path}/tick3.fastq.gz ${file_path}/tick4.fastq.gz ${file_path}/tick5.fastq.gz ${file_path}/tick6.fastq.gz ${file_path}/tick7.fastq.gz ${file_path}/tick8.fastq.gz > ${res_path}/logs/match_ku_mb_allticks.log

# Kraken 2
for t in tick1 tick2 tick3 tick4 tick5 tick6 tick7 tick8;
  do
    # Suppress a lot of output to be fair
    ./cgmemtime/cgmemtime ./k2/kraken2/k2 classify --threads 10 --db ./k2/standard_db --output - --report ./k2/tick_perf_k2_${t}.tsv ${file_path}/${t}.fastq.gz  > ${res_path}/logs/match_k2_std_${t}.log
  done

./cgmemtime/cgmemtime ./k2/kraken2/k2 classify --threads 10 --db ./k2/standard_db --output - --report ./k2/tick_perf_k2_allticks.tsv ${file_path}/tick1.fastq.gz ${file_path}/tick2.fastq.gz ${file_path}/tick3.fastq.gz ${file_path}/tick4.fastq.gz ${file_path}/tick5.fastq.gz ${file_path}/tick6.fastq.gz ${file_path}/tick7.fastq.gz ${file_path}/tick8.fastq.gz  > ${res_path}/logs/match_k2_std_allticks.log

done

# Ganon
for t in tick1 tick2 tick3 tick4 tick5 tick6 tick7 tick8;
  do
    # Suppress a lot of output to be fair
    ./cgmemtime/cgmemtime ganon classify --rel-cutoff 0 --rel-filter 0 --db-prefix ${basedir}/ganon/tick-borne_db -o ${basedir}/ganon/tick_perf_ganon_${t} --threads 10 -s ${file_path}/${t}.fastq.gz > ${res_path}/logs/match_ganon_viral_${t}.log
  done

./cgmemtime/cgmemtime ganon classify --rel-cutoff 0 --rel-filter 0 --db-prefix ${basedir}/ganon/tick-borne_db -o ${basedir}/ganon/tick_perf_ganon_allticks --threads 10 -s ${file_path}/tick2.fastq.gz ${file_path}/tick3.fastq.gz ${file_path}/tick4.fastq.gz ${file_path}/tick5.fastq.gz ${file_path}/tick6.fastq.gz ${file_path}/tick7.fastq.gz ${file_path}/tick8.fastq.gz > ${res_path}/logs/match_ganon_viral_allticks.log
