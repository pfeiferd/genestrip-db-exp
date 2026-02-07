#!/bin/sh
set -e

scriptdir=$(dirname "$0")

cd $scriptdir/..

mvn exec:exec@match -Dname=human_virus -Dgoal=extractrefseqfasta
mvn exec:exec@match -Dname=viral -Dgoal=extractrefseqfasta
