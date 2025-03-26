#!/bin/sh
set -e

scriptdir=$(dirname "$0")

cd $scriptdir/..

mkdir ku
cd ku

# We simply use the latest GitHub version... (currently v1.0.4)
git clone https://github.com/fbreitwieser/krakenuniq.git

cd krakenuniq

./install_krakenuniq -j .

# Create a viral database
mkdir ../viral_db
./krakenuniq-download --db ../viral_db refseq/viral
./krakenuniq-download --db ../viral_db taxonomy
./krakenuniq-build --db ../viral_db



