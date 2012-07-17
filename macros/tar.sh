#!/bin/sh

# Suppose you want to work with images of type 'TIF'...
#
# To compress image directory aaa/bbb/TIF and put result at ccc/
#
# > ./tar.csh -c TIF aaa/bbb/TIF ccc
#
# To uncompress image directory aaa/bbb/TIF and put result at ccc/
#
# > ./tar.csh -u TIF aaa/bbb/TIF ccc
#


if [ $(printf "%.1s" $4) = "/" ]
then
	dsttop=$4
else
	dsttop=$PWD/$4
fi

if [ "$1" = "-c" ]
then # compress

	dsttop=$dsttop/$2"_gzip"/$2

	cd $3

	for srcdir in $(ls -d V91*)
	do
		echo $srcdir
		dstdir=$dsttop/$srcdir
		mkdir -p $dstdir
		tar -zcf $dstdir/Plate1_0.tar.gz $srcdir/Plate1_0
	done

elif [ "$1" = "-u" ]
then # uncompress

	dsttop=$dsttop/$2"_unzp"/$2
	mkdir -p $dsttop

	cd $3

	for srcdir in $(ls -d V91*)
	do
		echo $srcdir
		tar -zxf $srcdir/Plate1_0.tar.gz -C $dsttop
	done
else
	echo "Error: Arg #1 must be -c or -u"
	exit
fi


#qsub -N del-$srcdir -cwd -V -b y -pe batch 8 "tar -zcf $dstdir/Plate1_0.tar.gz $srcdir/Plate1_0"

#qsub -N del-$srcdir -cwd -V -b y -pe batch 8 "tar -zxf $srcdir/Plate1_0.tar.gz -C $dsttop"


