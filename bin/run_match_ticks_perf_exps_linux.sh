#!/bin/sh
set -e

scriptdir=$(dirname "$0")

cd $scriptdir/..

for t in tick1 tick2 tick3 tick4 tick5 tick6 tick7 tick8;
  do
    rm -f ./data/projects/tick-borne/csv/tick-borne_match_${t}.csv
    ./cgmemtime/cgmemtime mvn exec:exec@singlematch -Dname=tick-borne -Dgoal=match -Dfqfile=${t}.fastq.gz  >& data/match_tick-borne_${t}.log
  done

rm -f ./data/projects/tick-borne/csv/tick-borne_match_allticks.csv
./cgmemtime/cgmemtime mvn exec:exec@match -Dname=tick-borne -Dgoal=match -Dfqmap=allticks.txt  >& data/match_tick-borne_allticks.log
