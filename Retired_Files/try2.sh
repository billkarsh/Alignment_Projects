#!/bin/sh


if [ $# = 1 ]
then
    last=$1
else
    last=$2
fi

for i in `seq $1 $last`
do
	j=$(($i + 1))
	echo $i $j
	./combine $i $j
	cd stack
	lsq pts.all ../../ldir -scale=1000 -square=1000 > log
	egrep 'Bad pair|RMS' < log
	cd ..
done


