#!/bin/sh
set -e

scriptdir=$(dirname "$0")

cd $scriptdir/..

for db in vineyard parasites human_virus viral tick-borne protozoa;
  do
    mvn exec:exec@db -Dname=$db -Dgoal=clear
    mvn exec:exec@db -Dname=$db -Dgoal=refseqfna
    mvn exec:exec@db -Dname=$db -Dgoal=fastasgenbankdl
    ./cgmemtime/cgmemtime -l mvn exec:exec@db -Dname=$db -Dgoal=db  >& data/db_gen_${db}.log
    mvn exec:exec@db -Dname=$db -Dgoal=dbinfo
  done