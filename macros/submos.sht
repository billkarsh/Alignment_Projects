#!/bin/sh

export MRC_TRIM=12

if (($# == 1))
then
	last=$1
else
	last=$2
fi

for lyr in $(seq $1 $last)
do
	echo $lyr
	if [ -d "$lyr" ]
	then
		qsub -N mos-$i -cwd -V -b y -pe batch 8 "mos simple 0,0,-1,-1 $i,$i -warp -nf > mos_$i.txt"
	fi
done

