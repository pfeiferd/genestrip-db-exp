#!/bin/sh
set -e

scriptdir=$(dirname "$0")

cd $scriptdir/..

mkdir -p sra-tools
cd sra-tools

wget https://github.com/gsauthof/cgmemtime/archive/refs/heads/master.zip
unzip master.zip
mv  cgmemtime-master  cgmemtime
cd  cgmemtime
make
