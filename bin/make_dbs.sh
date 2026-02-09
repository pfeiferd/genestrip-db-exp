#!/bin/sh
set -e

scriptdir=$(dirname "$0")

mkdir -p $scriptdir/..

cd $scriptdir/..

# Prepare genomes for other systems via Genestrip
mvn exec:exec@match -Dname=tick-borne -Dgoal=extractrefseqcsv
mvn exec:exec@match -Dname=human_virus -Dgoal=extractrefseqcsv
mvn exec:exec@match -Dname=viral -Dgoal=extractrefseqcsv

# Prepare input files for custom DB builds of other systems.
mvn exec:exec@inputcsv


### Ganon ###

# Create a dirs
mkdir -p ../ganon/viral_db
mkdir -p ../ganon/viral_lowfp_db
mkdir -p ../ganon/human_db
mkdir -p ../ganon/human_lowfp_db
mkdir -p ../ganon/tick-borne_db
mkdir -p ../ganon/tick-borne_lowfp_db

# Build ganon databases
ganon build-custom --input-file ../data/projects/viral/csv/viral_ganon.tsv --taxonomy-files ../data/common/nodes.dmp ../data/common/names.dmp --db-prefix ../ganon/viral_db --level species --threads 32
ganon build-custom --input-file ../data/projects/viral/csv/viral_ganon.tsv --taxonomy-files ../data/common/nodes.dmp ../data/common/names.dmp --db-prefix ../ganon/viral_lowfp_db --level species --threads 32 --max-fp 0.0000001
ganon build-custom --input-file ../data/projects/viral/csv/human_virus_ganon.tsv --taxonomy-files ../data/common/nodes.dmp ../data/common/names.dmp --db-prefix ../ganon/human_virus_db --level species --threads 32
ganon build-custom --input-file ../data/projects/viral/csv/human_virus_ganon.tsv --taxonomy-files ../data/common/nodes.dmp ../data/common/names.dmp --db-prefix ../ganon/human_virus_lowfp_db --level species --threads 32 --max-fp 0.0000001
ganon build-custom --input-file ../data/projects/viral/csv/tick-borne_ganon.tsv --taxonomy-files ../data/common/nodes.dmp ../data/common/names.dmp --db-prefix ../ganon/tick-borne_db --level species --threads 32
ganon build-custom --input-file ../data/projects/viral/csv/tick-borne_ganon.tsv --taxonomy-files ../data/common/nodes.dmp ../data/common/names.dmp --db-prefix ../ganon/tick-borne_lowfp_db --level species --threads 32 --max-fp 0.0000001

### KrakenUniq ###

cd $scriptdir/../ku/krakenuniq

# Create a viral database
mkdir -p ../viral_db
mkdir -p ../human_virus_db
mkdir -p ../tick-borne_db

mkdir -p ../viral_db/library
mkdir -p ../viral_db/taxonomy
mkdir -p ../human_virus_db/library
mkdir -p ../human_virus_db/taxonomy
mkdir -p ../tick-borne_db/library
mkdir -p ../tick-borne_db/taxonomy

cp ../../data/common/nodes.dmp ../viral_db/taxonomy
cp ../../data/common/names.dmp ../viral_db/taxonomy
cp ../../data/common/nodes.dmp ../human_virus_db/taxonomy
cp ../../data/common/names.dmp ../human_virus_db/taxonomy
cp ../../data/common/nodes.dmp ../tick-borne_db/taxonomy
cp ../../data/common/names.dmp ../tick-borne_db/taxonomy

cp ../../data/projects/viral/fasta/*.fa.gz ../viral_db/library
# KU only accepts unzipped files ending with fa, fna and so on
gunzip ../viral_db/library/*.fa.gz
cp ../../data/projects/viral/csv/viral_ku.map ../viral_db/library

cp ../../data/projects/human_virus/fasta/*.fa.gz ../human_virus_db/library
# KU only accepts unzipped files ending with fa, fna and so on
gunzip ../human_virus_db/library/*.fa.gz
cp ../../data/projects/human_virus_db/csv/human_virus_ku.map ../human_virus_db/library

cp ../../data/projects/tick-borne/fasta/*.fa.gz ../tick-borne_db/library
# KU only accepts unzipped files ending with fa, fna and so on
gunzip ../tick-borne_db/library/*.fa.gz
cp ../../data/projects/viral/csv/tick-borne_ku.map ../tick-borne_db/library

export JELLYFISH_BIN=$(pwd)/jellyfish-install/bin/jellyfish

../../cgmemtime/cgmemtime ./krakenuniq-build --db ../viral_db &> make_ku_viral_db_cgmemtime.log
../../cgmemtime/cgmemtime ./krakenuniq-build --db ../human_virus_db &> make_ku_human_virus_db_cgmemtime.log
../../cgmemtime/cgmemtime ./krakenuniq-build --db ../tick-borne_db &> make_ku_tick-borne_db_cgmemtime.log
