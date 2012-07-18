#!/bin/sh

# > ./subdel dirname

qsub -N del-$1 -cwd -V -b y -pe batch 8 rm -rf $1


