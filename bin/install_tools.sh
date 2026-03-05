#!/bin/sh
set -e

scriptdir=$(dirname "$0")

# iss and ganon are supposed to be installed manually via conda
# iss and ganon are supposed to be executable from the command line.
# For iss see:
# https://insilicoseq.readthedocs.io/en/latest/iss/install.html
# For ganon see:
# https://pirovc.github.io/ganon/start/#install

# gm - prerequisites on Ubuntu
echo "Preparing system... installing prerequisites"
sudo apt update
sudo apt install libbz2-dev zlib1g-dev

echo "DONE"
echo ""

### cgmemtime ###

cd $scriptdir/..
basedir=$(pwd)

# gm fix - only build cgmemtime if not already there
if [ ! -d "cgmemtime" ]; then
  wget https://github.com/gsauthof/cgmemtime/archive/refs/heads/master.zip
  unzip master.zip
  mv  cgmemtime-master  cgmemtime
  cd  cgmemtime
  # Let's print to stdout:
  sed -i 's/print_result(stderr/print_result(stdout/' cgmemtime.c
  make
  cd $scriptdir/..
  rm master.zip
fi

### sra-tools ###

cd ${basedir}

rm -rf sra-tools
mkdir -p sra-tools
cd sra-tools

wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/3.2.1/sratoolkit.3.2.1-alma_linux64.tar.gz
gunzip sratoolkit.3.2.1-alma_linux64.tar.gz
tar -xf sratoolkit.3.2.1-alma_linux64.tar
mv sratoolkit.3.2.1-alma_linux64 sratoolkit
rm sratoolkit.3.2.1-alma_linux64.tar

### kraken-uniq ###

cd ${basedir}

mkdir -p ku
cd ku

# We simply use the latest GitHub version... (currently v1.0.4)
rm -rf krakenuniq
git clone https://github.com/fbreitwieser/krakenuniq.git
# if error uint32_t appears, create a fixed version of uid_mapping.hpp with #include <cstdint> and report-cols.hpp in $basedir/bin/fixes/ku/ folder. 
# cp $basedir/bin/fixes/ku/*.hpp krakenuniq/src
cd krakenuniq

# This file must be deleted, it leads to a compile error as it is not C++
rm src/gzstream/version
# Fix missing include:
sed -i '1i\#include <stdint.h>' src/uid_mapping.hpp
sed -i '13s/class//' src/report-cols.hpp
# Now we can build...
./install_krakenuniq.sh -j .


### kraken 2 ###

cd ${basedir}

mkdir -p k2
cd k2

# We simply use the latest GitHub version... (currently v1.0.4)
git clone https://github.com/DerrickWood/kraken2.git

cd kraken2

# Now we can build...
./install_kraken2.sh .
