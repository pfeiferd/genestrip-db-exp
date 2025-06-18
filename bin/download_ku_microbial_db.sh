#!/bin/sh
set -e

scriptdir=$(dirname "$0")

cd $scriptdir/../ku/krakenuniq

# Create a viral database
mkdir -p ../microbial_db

cd ../microbial_db
wget https://genome-idx.s3.amazonaws.com/kraken/uniq/krakendb-2023-08-08-MICROBIAL/kuniq_microbialdb_minus_kdb.20230808.tgz
gunzip kuniq_microbialdb_minus_kdb.20230808.tgz
tar -xf kuniq_microbialdb_minus_kdb.20230808.tar
rm kuniq_microbialdb_minus_kdb.20230808.tar
wget https://genome-idx.s3.amazonaws.com/kraken/uniq/krakendb-2023-08-08-MICROBIAL/database.kdb

