Set Up
-------
Download latest version for sequential machines (here, superlu_4.3.tar.gz).

(In linux NX, just double-click gz file and drag folder to extract).
WARNING: If use Windows tools to extract the line endings may change to <CR><LF>.

Apply Bill_Changes/Sequential as follows...

Read SuperLU_4.3 ReadMe to edit make.inc, and to build libs.

*** IMPORTANT: Make sure editor uses unix line endings.

I made these changes to make.inc:
SuperLUroot		= /groups/apig/tomo/Libraries/SUPERLU/SuperLU_4.3
SUPERLULIB   	= $(SuperLUroot)/lib/libsuperlu_4.3.a
BLASLIB			= ../lib/blas$(PLAT).a

> make blaslib

######################################################
Note: Later I realized the cluster already has a BLAS library and I tried these lines in make.inc, and recompiled everything.

BLASDEF 	= -DUSE_VENDOR_BLAS
BLASLIB		= /usr/lib64/libblas.a

The result was identical.
######################################################

> make

Copy Bill_Changes/Sequential parts to appropriate subfolders.

In EXAMPLE/ directory:
> ./update.sht




==========================================================================
Set Up SuperLU_DIST_3.3
------------------------
Download latest version and using NX, drag out the SuperLU_DIST_3.3 folder.

Apply Bill_Changes/Distributed as follows...

Read README, and edit make.inc as follows:

Select MAKE_INC/make.i386_linux as base make.inc...

DSuperLUroot 	= /groups/apig/tomo/Libraries/SUPERLU/SuperLU_DIST_3.3
DSUPERLULIB   	= $(DSuperLUroot)/lib/libsuperlu_dist_3.3.a

BLASDEF	     	= -DUSE_VENDOR_BLAS
BLASLIB      	= /usr/lib64/libblas.a

Need METISLIB and PARMETISLIB, so follow README, get:
http://glaros.dtc.umn.edu/gkhome/metis/parmetis/download

Need CMake, so edit .bash_profile to include:
PATH=$PATH:/usr/local/cmake-2.8.8/bin

In parmetis-4.0.3 top directory:
> make config
> make

The make usually fails because command line option -Wno-unused-but-set-variable is not recognized by compiler. So go looking for all files 'flags.make' in the build directory and edit them to remove the flag.

Finally it will make, but libs cannot really be copied to /usr/local, so we edit the superLU make.inc thus:
METISLIB        = /groups/apig/tomo/Libraries/METIS/parmetis-4.0.3/build/Linux-x86_64/libmetis/libmetis.a
PARMETISLIB     = /groups/apig/tomo/Libraries/METIS/parmetis-4.0.3/build/Linux-x86_64/libparmetis/libparmetis.a

Now copy our files from Distributed/ to approriate folders.

> make

In EXAMPLE/ directory:
> ./update.sht







