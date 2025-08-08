#!/bin/sh
set -e

scriptdir=$(dirname "$0")

cd $scriptdir/..

wget https://github.com/gsauthof/cgmemtime/archive/refs/heads/master.zip
unzip master.zip
mv  cgmemtime-master  cgmemtime
cd  cgmemtime
# Let' print to stdout:
sed -i 's/print_result(stderr, &args, &res);/print_result(stdout, &args, &res);/g' cgmemtime.c
make

cd $scriptdir/..
rm master.zip
