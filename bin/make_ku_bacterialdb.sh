#!/bin/sh
set -e

scriptdir=$(dirname "$0")

cd $scriptdir/../ku/krakenuniq

# Create a viral database
mkdir ../bacteria_db
./krakenuniq-download --db ../bacteria_db --threads 10 refseq/bacteria
./krakenuniq-download --db ../bacteria_db taxonomy

export JELLYFISH_BIN=$scriptdir/../ku/krakenuniq/jellyfish-install/bin/jellyfish

./krakenuniq-build --db ../bacteria_db --threads 10 --work-on-disk
