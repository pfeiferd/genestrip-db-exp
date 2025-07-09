#!/bin/sh
set -e

scriptdir=$(dirname "$0")

cd $scriptdir/..

mvn exec:exec@match -Dname=human_virus -Dgoal=match -Dfqmap=saliva.txt
