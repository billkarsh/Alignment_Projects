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
		qsub -N lou-s-$lyr -cwd -V -b y -pe batch 8 make -f make.same -j 8 EXTRA='""'
		if (($lyr > $1))
		then
			qsub -N lou-d-$lyr -cwd -V -b y -pe batch 8 make -f make.down -j 8 EXTRA='""'
		fi
		cd ..
	fi
done

