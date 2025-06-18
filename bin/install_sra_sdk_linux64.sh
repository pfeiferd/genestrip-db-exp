#!/bin/sh
set -e

scriptdir=$(dirname "$0")

cd $scriptdir/..

mkdir -p sra-tools
cd sra-tools

wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/3.2.1/sratoolkit.3.2.1-alma_linux64.tar.gz
gunzip sratoolkit.3.2.1-alma_linux64.tar.gz
tar -xf sratoolkit.3.2.1-alma_linux64.tar.gz
mv sratoolkit.3.2.1-alma_linux64 sratoolkit
rm sratoolkit.3.2.1-alma_linux64.tar