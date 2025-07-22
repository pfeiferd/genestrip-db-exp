#!/bin/sh
set -e

scriptdir=$(dirname "$0")

cd $scriptdir/../ku/krakenuniq

# Create a viral database
mkdir -p ../viral_db
./krakenuniq-download --db ../viral_db refseq/viral/Any
./krakenuniq-download --db ../viral_db taxonomy
# Remove this single genome as one of its k-mers messes up the validation:
#mkdir ../viral_db/library/viral/moved
#mv ../viral_db/library/viral/Chromosome/Tenuivirus_zeae_na-tax3052767* ../viral_db/library/viral/moved
#sed '/Tenuivirus_zeae_na-tax3052767/d' ../viral_db/library-files.txt > ../viral_db/library-files.txt.new
#mv ../viral_db/library-files.txt ../viral_db/library-files.txt.org
#mv ../viral_db/library-files.txt.new ../viral_db/library-files.txt

export JELLYFISH_BIN=$(pwd)/jellyfish-install/bin/jellyfish

./krakenuniq-build --db ../viral_db
