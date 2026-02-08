#!/bin/sh
set -e

scriptdir=$(dirname "$0")

cd $scriptdir/..

mvn exec:exec@match -Dname=tick-borne -Dgoal=extractrefseqcsv
mvn exec:exec@match -Dname=human_virus -Dgoal=extractrefseqcsv
mvn exec:exec@match -Dname=viral -Dgoal=extractrefseqcsv
