#!/bin/bash
set -e

scriptdir=$(dirname "$0")

cd $scriptdir/..

res_path=./results

for t in tick1 tick2 tick3 tick4 tick5 tick6 tick7 tick8;
  do
    rm -f ./data/projects/tick-borne/csv/tick-borne_match_${t}.csv
    ./cgmemtime/cgmemtime mvn exec:exec@singlematch -Dname=tick-borne -Dgoal=match -Dfqfile=${t}.fastq.gz  >& ${res_path}/logs/match_tick-borne_${t}.log
    mv ./data/projects/tick-borne/csv/tick-borne_match_${t}.csv ${res_path}/reports
  done

./cgmemtime/cgmemtime mvn exec:exec@match -Dname=tick-borne -Dgoal=match -Dfqmap=ticks.txt  >& ${res_path}/logs/match_tick-borne_allticks.log
