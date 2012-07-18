#!/bin/sh

rm -f pts.all

#get line 1, subst 'IDBPATH=xxx' with 'xxx'
idb=$(sed -n -e 's|IDBPATH \(.*\)|\1|' -e '1p' < imageparams.txt)

cp imageparams.txt pts.all

for i in $(seq $1 $2)
do
	cat $idb/$i/fm.same >> pts.all
done

for i in $(seq $1 $2)
do
	echo $i
	if (($i == $1))
	then
		cat $i/pts.{same} >> pts.all
	else
		cat $i/pts.{down,same} >> pts.all
	fi
done

mv pts.all stack

