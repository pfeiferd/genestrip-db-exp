#!/bin/sh
set -e

scriptdir=$(dirname "$0")

cd $scriptdir/..
basedir=$(pwd)

# NONE of the following is needed (yet) - it does not get execute since the x-loop is empty.
for x in ; do
for db in standard;
  # Kraken 2
  do
   for id in tick1_sim tick2_sim tick3_sim tick4_sim tick5_sim tick6_sim tick7_sim tick8_sim;
    do
      ./k2/kraken2/k2 classify --threads 10 --db ./k2/${db}_db ./data/fastq/${id}.fastq --output ./k2/${db}_${id}.tsv
      # Same again with high confidence of 0.8
      ./k2/kraken2/k2 classify --confidence 0.8 --threads 10 --db ./k2/${db}_db ./data/fastq/${id}.fastq --output ./k2/${db}_highconf_${id}.tsv
    done
  done

  ## Genestrip ##
  mvn exec:exec@matchrep -Dname=tick-borne -Dgoal=match -Dfqmap=ticks.txt
  mvn exec:exec@matchrep -Dname=tick-borne -Dgoal=match -Dfqmap=ticks_sim.txt
done


## Ganon ##

# Make sure Nanosim is deactivated as it is associated with an older version of Ganon!
# conda deactivate

for db in tick-borne; # tick-borne_lowfp;
  do

  # Ganon on simulated tick files
  for id in tick1_sim tick2_sim tick3_sim tick4_sim tick5_sim tick6_sim tick7_sim tick8_sim;
    do
      ganon classify --db-prefix ${basedir}/ganon/${db}_db -s ./data/fastq/${id}.fastq --output-all -o ./ganon/${db}_${id} --threads 32 --fpr-query 1e-0 --rel-cutoff 0 --rel-filter 0
    done

  # Ganon on tick files
#  for id in tick1 tick2 tick3 tick4 tick5 tick6 tick7 tick8;
#    do
#      ganon classify --db-prefix ./ganon/${db}_db -s ./data/fastq/${id}.fastq.gz --output-all -o ./ganon/${db}_${id} --threads 32
#    done
  done

