#!/bin/sh
set -e

scriptdir=$(dirname "$0")

cd $scriptdir/..

mvn exec:exec@acc
mvn exec:exec@exp
