#!/bin/sh
set -e

scriptdir=$(dirname "$0")

cd $scriptdir/..

/usr/bin/time -l mvn exec:exec@exp
