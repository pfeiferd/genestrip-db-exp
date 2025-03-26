#!/bin/sh
set -e

scriptdir=$(dirname "$0")

cd $scriptdir/..

mkdir -p ku
cd ku

# We simply use the latest GitHub version... (currently v1.0.4)
git clone https://github.com/fbreitwieser/krakenuniq.git

cd krakenuniq

./install_krakenuniq.sh -j .



