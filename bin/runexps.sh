#!/bin/sh
set -e

scriptdir=$(dirname "$0")

cd $scriptdir/..

mvn exec:exec@exp
