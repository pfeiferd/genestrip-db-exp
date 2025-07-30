#!/bin/sh
set -e

scriptdir=$(dirname "$0")

cd $scriptdir/..

wget https://github.com/gsauthof/cgmemtime/archive/refs/heads/master.zip
unzip master.zip
mv  cgmemtime-master  cgmemtime
cd  cgmemtime
make

cd $scriptdir/..
rm master.zip
