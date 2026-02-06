#!/bin/sh
set -e

scriptdir=$(dirname "$0")

cd $scriptdir/..

mkdir -p ku
cd ku

# We simply use the latest GitHub version... (currently v1.0.4)
git clone https://github.com/fbreitwieser/krakenuniq.git

cd krakenuniq

# This file must be deleted, it leads to a compile error as it is not C++
rm src/gzstream/version
# Now we can build...
./install_krakenuniq.sh -j .



