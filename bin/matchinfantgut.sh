#!/bin/sh
set -e

scriptdir=$(dirname "$0")

cd $scriptdir/..

if [ -z "$1" ]
  then
	goal=match
  else
  	goal=$1
fi

mvn exec:exec@match -Dfqmap=./data/fastq/infant_gut.txt -Dname=human_gut -Dgoal=$goal
mvn exec:exec@match -Dfqmap=./data/fastq/infant_gut.txt -Dname=human_gut_fungi -Dgoal=$goal
