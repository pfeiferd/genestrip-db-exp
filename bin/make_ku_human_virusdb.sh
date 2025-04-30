#!/bin/sh
set -e

scriptdir=$(dirname "$0")

cd $scriptdir/../ku/krakenuniq

# Create a human_virus database
mkdir -p ../human_virus_db

file_path=../../data/projects/human_virus/unfolded_taxids.txt

line_number=1
while read -r line; do
      ./krakenuniq-download --db ../human_virus_db refseq/viral/Any/taxid=$line
done < "$file_path"

./krakenuniq-download --db ../human_virus_db taxonomy

export JELLYFISH_BIN=$(pwd)/jellyfish-install/bin/jellyfish

./krakenuniq-build --db ../human_virus_db
