#!/bin/sh
set -e

scriptdir=$(dirname "$0")

cd $scriptdir/..
basedir=$(pwd)

if [0 == 1]; then
# Prepare genomes for other systems via Genestrip
mvn exec:exec@match -Dname=human_virus -Dgoal=extractrefseqcsv
mvn exec:exec@match -Dname=viral -Dgoal=extractrefseqcsv
mvn exec:exec@match -Dname=tick-borne -Dgoal=extractrefseqcsv

# Prepare input files for custom DB builds of other systems.
mvn exec:exec@inputcsv

#### Ganon ###

# Create a dirs
mkdir -p ganon/viral_db
mkdir -p ganon/viral_lowfp_db
mkdir -p ganon/human_virus_db
mkdir -p ganon/human_lowfp_db
# Not needed (yet):
#mkdir -p ganon/tick-borne_db
#mkdir -p ganon/tick-borne_lowfp_db

# Build ganon databases
ganon build-custom --input-file data/projects/viral/csv/viral_ganon.tsv --taxonomy-files data/common/nodes.dmp data/common/names.dmp --db-prefix ganon/viral_db --level leaves --threads 32
ganon build-custom --input-file data/projects/viral/csv/viral_ganon.tsv --taxonomy-files data/common/nodes.dmp data/common/names.dmp --db-prefix ganon/viral_lowfp_db --level leaves --threads 32 --max-fp 0.0000001
ganon build-custom --input-file data/projects/human_virus/csv/human_virus_ganon.tsv --taxonomy-files data/common/nodes.dmp data/common/names.dmp --db-prefix ganon/human_virus_db --level leaves --threads 32
ganon build-custom --input-file data/projects/human_virus/csv/human_virus_ganon.tsv --taxonomy-files data/common/nodes.dmp data/common/names.dmp --db-prefix ganon/human_virus_lowfp_db --level leaves --threads 32 --max-fp 0.0000001
## Not needed (yet):
##ganon build-custom --input-file data/projects/tick-borne/csv/tick-borne_ganon.tsv --taxonomy-files data/common/nodes.dmp data/common/names.dmp --db-prefix ganon/tick-borne_db --level leaves --threads 32
##ganon build-custom --input-file data/projects/tick-borne/csv/tick-borne_ganon.tsv --taxonomy-files data/common/nodes.dmp data/common/names.dmp --db-prefix ganon/tick-borne_lowfp_db --level leaves --threads 32 --max-fp 0.0000001

fi
### KrakenUniq ###

cd $basedir/ku/krakenuniq

# Create a viral database
mkdir -p ../viral_db
mkdir -p ../human_virus_db
mkdir -p ../tick-borne_db

mkdir -p ../viral_db/library
mkdir -p ../viral_db/taxonomy
mkdir -p ../human_virus_db/library
mkdir -p ../human_virus_db/taxonomy
# Not needed (yet):
#mkdir -p ../tick-borne_db/library
#mkdir -p ../tick-borne_db/taxonomy

cp ../../data/common/nodes.dmp ../viral_db/taxonomy
cp ../../data/common/names.dmp ../viral_db/taxonomy
cp ../../data/common/nodes.dmp ../human_virus_db/taxonomy
cp ../../data/common/names.dmp ../human_virus_db/taxonomy
# Not needed (yet):
#cp ../../data/common/nodes.dmp ../tick-borne_db/taxonomy
#cp ../../data/common/names.dmp ../tick-borne_db/taxonomy

# KU only accepts unzipped files ending with fa, fna and so on
cp ../../data/projects/viral/fasta/*.fa ../viral_db/library
cp ../../data/projects/viral/csv/viral_ku.map ../viral_db/library

cp ../../data/projects/human_virus/fasta/*.fa ../human_virus_db/library
cp ../../data/projects/human_virus/csv/human_virus_ku.map ../human_virus_db/library

# Not needed (yet):
#cp ../../data/projects/tick-borne/fasta/*.fa ../tick-borne_db/library
#cp ../../data/projects/viral/csv/tick-borne_ku.map ../tick-borne_db/library

export JELLYFISH_BIN=$(pwd)/jellyfish-install/bin/jellyfish

./krakenuniq-build --db ../viral_db # &> make_ku_viral_db_cgmemtime.log
./krakenuniq-build --db ../human_virus_db # &> make_ku_human_virus_db_cgmemtime.log
# Not needed (yet):
# ../../cgmemtime/cgmemtime ./krakenuniq-build --db ../tick-borne_db # &> make_ku_tick-borne_db_cgmemtime.log

### Kraken 2 ###

cd $basedir/k2/kraken2

# Create a viral database
mkdir -p ../viral_db
mkdir -p ../human_virus_db
mkdir -p ../tick-borne_db

mkdir -p ../viral_db/library
mkdir -p ../viral_db/taxonomy
mkdir -p ../human_virus_db/library
mkdir -p ../human_virus_db/taxonomy
# Not needed (yet):
#mkdir -p ../tick-borne_db/library
#mkdir -p ../tick-borne_db/taxonomy

cp ../../data/common/nodes.dmp ../viral_db/taxonomy
cp ../../data/common/names.dmp ../viral_db/taxonomy
cp ../../data/common/nodes.dmp ../human_virus_db/taxonomy
cp ../../data/common/names.dmp ../human_virus_db/taxonomy
# Not needed (yet):
#cp ../../data/common/nodes.dmp ../tick-borne_db/taxonomy
#cp ../../data/common/names.dmp ../tick-borne_db/taxonomy

for file in $basedir/ku/viral_db/library/*.fa
do
    ./kraken2-build --no-masking --add-to-library $file --db ../viral_db
done
./kraken2-build --no-masking --build --db ../viral_db

for file in $basedir/ku/human_virus_db/library/*.fa
do
    ./kraken2-build --no-masking --add-to-library $file --db ../human_virus_db
done
./kraken2-build --no-masking --build --db ../human_virus_db


## Download ready made databases for tick analysis:

### KrakenUniq microbial-db ###

cd $basedir/ku/krakenuniq

mkdir -p ../microbial_db

cd ../microbial_db
wget https://genome-idx.s3.amazonaws.com/kraken/uniq/krakendb-2023-08-08-MICROBIAL/kuniq_microbialdb_minus_kdb.20230808.tgz
gunzip kuniq_microbialdb_minus_kdb.20230808.tgz
tar -xf kuniq_microbialdb_minus_kdb.20230808.tar
rm kuniq_microbialdb_minus_kdb.20230808.tar
wget https://genome-idx.s3.amazonaws.com/kraken/uniq/krakendb-2023-08-08-MICROBIAL/database.kdb

### Kraken 2

cd $basedir/k2/kraken2

mkdir -p ../standard_db

cd ../standard_db
wget https://genome-idx.s3.amazonaws.com/kraken/k2_standard_20251015.tar.gz
gunzip k2_standard_20251015.tar.gz
tar -xf k2_standard_20251015.tar
rm k2_standard_20251015.tar

### Ganon

cd $basedir

mkdir -p ganon/standard_db

ganon build --source refseq --organism-group bacteria --threads 48 --complete-genomes --db-prefix ganon/standard_db