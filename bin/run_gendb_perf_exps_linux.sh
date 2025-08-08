#!/bin/sh
set -e

scriptdir=$(dirname "$0")

cd $scriptdir/..

res_path=./results

for db in vineyard parasites human_virus viral tick-borne protozoa;
  do
    mvn exec:exec@db -Dname=$db -Dgoal=clear
#    mvn exec:exec@db -Dname=$db -Dgoal=refseqfna
#    mvn exec:exec@db -Dname=$db -Dgoal=fastasgenbankdl
    ./cgmemtime/cgmemtime mvn exec:exec@db -Dname=$db -Dgoal=db  > ${res_path}/logs/db_gen_${db}.log
    mvn exec:exec@db -Dname=$db -Dgoal=dbinfo
  done