#!/bin/csh

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


if( "$4" =~ "/*" ) then
	set dsttop = $4
else
	set dsttop = $PWD/$4
endif


if( "$1" =~ "-c" ) then # compress

	set dsttop = $dsttop/$2"_gzip"/$2

	cd $3

	foreach srcdir (`ls -d V91*`)
		echo $srcdir
		set dstdir = $dsttop/$srcdir
		mkdir -p $dstdir
		tar -zcf $dstdir/Plate1_0.tar.gz $srcdir/Plate1_0
	end

else if( "$1" =~ "-u" ) then # uncompress

	set dsttop = $dsttop/$2"_unzp"/$2
	mkdir -p $dsttop

	cd $3

	foreach srcdir (`ls -d V91*`)
		echo $srcdir
		tar -zxf $srcdir/Plate1_0.tar.gz -C $dsttop
	end
else
	echo "Error: Arg #1 must be -c or -u"
	exit
endif


#qsub -N del-$srcdir -cwd -V -b y -pe batch 8 "tar -zcf $dstdir/Plate1_0.tar.gz $srcdir/Plate1_0"

#qsub -N del-$srcdir -cwd -V -b y -pe batch 8 "tar -zxf $srcdir/Plate1_0.tar.gz -C $dsttop"


