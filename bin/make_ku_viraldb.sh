#!/bin/sh
set -e

scriptdir=$(dirname "$0")

cd $scriptdir/../ku/krakenuniq

# Create a viral database
mkdir -p ../viral_db
# ./krakenuniq-download --db ../viral_db refseq/viral/Any
#./krakenuniq-download --db ../viral_db taxonomy

mkdir -p ../viral_db/library
mkdir -p ../viral_db/taxonomy

cp ../../data/common/nodes.dmp ../viral_db/taxonomy
cp ../../data/common/names.dmp ../viral_db/taxonomy

cp ../../data/projects/viral/fasta/*.fa.gz ../viral_db/library
# KU only accepts unzipped files ending with fa, fna and so on
gunzip ../../data/projects/viral/fasta/*.fa.gz
cp ../../data/projects/viral/csv/viral_ku.map ../viral_db/library

export JELLYFISH_BIN=$(pwd)/jellyfish-install/bin/jellyfish

../../cgmemtime/cgmemtime ./krakenuniq-build --db ../viral_db &> make_ku_viraldb_cgmemtime.log
