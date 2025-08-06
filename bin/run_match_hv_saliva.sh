#!/bin/sh
set -e

scriptdir=$(dirname "$0")

cd $scriptdir/..

res_path=./results

mvn exec:exec@match -Dname=human_virus -Dgoal=match -Dfqmap=saliva.txt

mv ./data/projects/viral/csv/human_virus_*.csv ${res_path}/reports
