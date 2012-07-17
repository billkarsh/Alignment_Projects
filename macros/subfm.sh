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
		cd $lyr
		qsub -N lou-f-$lyr -cwd -V -b y -pe batch 8 make -f make.fm -j 8 EXTRA='""'
		cd ..
	fi
done

