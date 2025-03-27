#!/bin/sh
set -e

scriptdir=$(dirname "$0")

cd $scriptdir/../ku/krakenuniq

# Create a viral database
mkdir ../viral_db
./krakenuniq-download --db ../viral_db refseq/viral
./krakenuniq-download --db ../viral_db taxonomy

export JELLYFISH_BIN=$scriptdir/../ku/krakenuniq/jellyfish-install/bin

./krakenuniq-build --db ../viral_db
