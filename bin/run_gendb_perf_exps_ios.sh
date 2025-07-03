#!/bin/sh
set -e

scriptdir=$(dirname "$0")

cd $scriptdir/..

if [ -z "$1" ]
  then
	goal=db
  else
  	goal=$1
fi

/usr/bin/time -l mvn exec:exec@db -Dname=human_virus -Dgoal=$goal  >& data/db_gen_hv.log
mvn exec:exec@db -Dname=human_virus -Dgoal=dbinfo
#/usr/bin/time -l mvn exec:exec@db -Dname=virus -Dgoal=$goal  >& db_gen_cv.log
#/usr/bin/time -l mvn exec:exec@db -Dname=tick-borne -Dgoal=$goal  >& db_gen_tb.log
#/usr/bin/time -l mvn exec:exec@db -Dname=vineyard -Dgoal=$goal  >& db_gen_vin.log
#/usr/bin/time -l mvn exec:exec@db -Dname=protozoa -Dgoal=$goal  >& db_gen_pro.log
#/usr/bin/time -l mvn exec:exec@db -Dname=parasites -Dgoal=$goal  >& db_gen_par.log
