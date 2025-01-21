#!/bin/sh
set -e

scriptdir=$(dirname "$0")

cd $scriptdir/..

if [ -z "$1" ]
  then
	goal=genall
  else
  	goal=$1
fi


mvn exec:exec@db -Dname=human_virus2 -Dgoal=$goal
mvn exec:exec@db -Dname=human_virus2-minupdate -Dgoal=$goal
mvn exec:exec@db -Dname=babesia -Dgoal=$goal
mvn exec:exec@db -Dname=babesia-minupdate -Dgoal=$goal
mvn exec:exec@db -Dname=protozoa -Dgoal=$goal
mvn exec:exec@db -Dname=protozoa-minupdate -Dgoal=$goal
mvn exec:exec@db -Dname=protozoa-rna -Dgoal=$goal
mvn exec:exec@db -Dname=protozoa-rna-minupdate -Dgoal=$goal
mvn exec:exec@db -Dname=vineyard -Dgoal=$goal
mvn exec:exec@db -Dname=vineyard-minupdate -Dgoal=$goal
mvn exec:exec@db -Dname=parasites -Dgoal=$goal
mvn exec:exec@db -Dname=parasites-minupdate -Dgoal=$goal
mvn exec:exec@db -Dname=borrelia -Dgoal=$goal
mvn exec:exec@db -Dname=borrelia-minupdate -Dgoal=$goal
mvn exec:exec@db -Dname=borrelia_plasmid -Dgoal=$goal
mvn exec:exec@db -Dname=borrelia_plasmid-minupdate -Dgoal=$goal
mvn exec:exec@db -Dname=chronicb -Dgoal=$goal
mvn exec:exec@db -Dname=chronicb-minupdate -Dgoal=$goal
mvn exec:exec@db -Dname=chronicb-rna -Dgoal=$goal
mvn exec:exec@db -Dname=chronicb-rna-minupdate -Dgoal=$goal
mvn exec:exec@db -Dname=plasmopara -Dgoal=$goal
mvn exec:exec@db -Dname=plasmopara-minupdate -Dgoal=$goal
mvn exec:exec@db -Dname=human_gut -Dgoal=$goal

