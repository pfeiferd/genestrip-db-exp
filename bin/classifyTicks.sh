#!/bin/sh
set -e

scriptdir=$(dirname "$0")

cd $scriptdir/..
basedir=$(pwd)


# General note: We skipped the *_2.fastq.gz files from the human saliva analysis because the files are so huge
# and analysis takes too long - in particular under ganon...

## Genestrip ##

mvn exec:exec@matchrep2 -Dname=tick-borne -Dgoal=match -Dfqmap=ticks.txt
mvn exec:exec@matchrep2 -Dname=tick-borne -Dgoal=match -Dfqmap=ticks_sim.txt

## Ganon ##

for x in ; do

# Ganon on simulated tick files
for id in tick1_sim #tick2 tick2_sim ...
  do
    ganon classify --db-prefix ./ganon/standard_db -s ./data/fastq/${id}.fastq --output-all -o ./ganon/${id} --threads 32
  done

# Ganon on tick files
for id in tick1 #tick2 ...
  do
    ganon classify --db-prefix ./ganon/standard_db -s ./data/fastq/${id}.fastq.gz --output-all -o ./ganon/${id} --threads 32
  done

done
