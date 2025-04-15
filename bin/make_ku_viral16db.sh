#!/bin/sh
set -e

scriptdir=$(dirname "$0")

cd $scriptdir/../ku/krakenuniq

# Create a viral database
mkdir -p ../viral_db
./krakenuniq-download --db ../viral_db refseq/viral/Any
./krakenuniq-download --db ../viral_db taxonomy

export JELLYFISH_BIN=$(pwd)/jellyfish-install/bin/jellyfish

./krakenuniq-build --db ../viral_db
