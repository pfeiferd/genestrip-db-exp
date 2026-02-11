#!/bin/sh
set -e

scriptdir=$(dirname "$0")

cd $scriptdir/..
basedir=$(pwd)

cd $basedir/k2/kraken2

for file in $basedir/ku/viral_db/library/*.fa
do
    ./kraken2-build --add-to-library $file --db ../viral_db
done
./kraken2-build --build --db ../viral_db

for file in $basedir/ku/human_virus_db/library/*.fa
do
    ./kraken2-build --add-to-library $file --db ../human_virus_db
done
./kraken2-build --build --db ../human_virus_db
