#!/bin/sh
set -e

scriptdir=$(dirname "$0")

cd $scriptdir/../ku/krakenuniq

# Create a viral database
mkdir -p ../viral16_db
./krakenuniq-download --db ../viral16_db refseq/viral/Any
./krakenuniq-download --db ../viral16_db taxonomy

export JELLYFISH_BIN=$(pwd)/jellyfish-install/bin/jellyfish

./krakenuniq-build --kmer-len 16 --db ../viral16_db

