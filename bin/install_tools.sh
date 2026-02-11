#!/bin/sh
set -e

scriptdir=$(dirname "$0")

# iss and ganon are supposed to be installed manually via conda
# iss and ganon are supposed to be executable from the command line.
# For iss see:
# https://insilicoseq.readthedocs.io/en/latest/iss/install.html
# For ganon see:
# https://pirovc.github.io/ganon/start/#install

### cgmemtime ###

cd $scriptdir/..

wget https://github.com/gsauthof/cgmemtime/archive/refs/heads/master.zip
unzip master.zip
mv  cgmemtime-master  cgmemtime
cd  cgmemtime
# Let's print to stdout:
sed -i 's/print_result(stderr/print_result(stdout/' cgmemtime.c
make

cd $scriptdir/..
rm master.zip

### sra-tools ###

cd $scriptdir/..

mkdir -p sra-tools
cd sra-tools

wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/3.2.1/sratoolkit.3.2.1-alma_linux64.tar.gz
gunzip sratoolkit.3.2.1-alma_linux64.tar.gz
tar -xf sratoolkit.3.2.1-alma_linux64.tar
mv sratoolkit.3.2.1-alma_linux64 sratoolkit
rm sratoolkit.3.2.1-alma_linux64.tar

### kraken-uniq ###

mkdir -p ku
cd ku

# We simply use the latest GitHub version... (currently v1.0.4)
git clone https://github.com/fbreitwieser/krakenuniq.git

cd krakenuniq

# This file must be deleted, it leads to a compile error as it is not C++
rm src/gzstream/version
# Now we can build...
./install_krakenuniq.sh -j .


### kraken 2 ###

mkdir -p k2
cd k2

# We simply use the latest GitHub version... (currently v1.0.4)
git clone https://github.com/DerrickWood/kraken2.git

cd kraken2

# Now we can build...
./install_kraken2.sh .
