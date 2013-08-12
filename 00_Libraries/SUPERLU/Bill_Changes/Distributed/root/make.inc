############################################################################
#
#  Program:         SuperLU_DIST
#
#  Module:          make.inc
#
#  Purpose:         Top-level Definitions
#
#  Creation date:   February 4, 1999   version alpha
#
#  Modified:	    September 1, 1999  version 1.0
#                   March 15, 2003     version 2.0
#
# 		    January 18, 2006   Sam Adams
#                                      General Dynamics - Network Systems
#                       works for i386 Linux, with LAM-MPI 7.1.1 and GCC 4.
#
############################################################################
#
#  The machine (platform) identifier to append to the library names
#
PLAT		= _i386

#
#  The name of the libraries to be created/linked to
#
DSuperLUroot 	= /groups/apig/tomo/Libraries/SUPERLU/SuperLU_DIST_3.3
DSUPERLULIB   	= $(DSuperLUroot)/lib/libsuperlu_dist_3.3.a
#
BLASDEF	     	= -DUSE_VENDOR_BLAS
BLASLIB      	= /usr/lib64/libblas.a
METISLIB        = /groups/apig/tomo/Libraries/METIS/parmetis-4.0.3/build/Linux-x86_64/libmetis/libmetis.a
PARMETISLIB     = /groups/apig/tomo/Libraries/METIS/parmetis-4.0.3/build/Linux-x86_64/libparmetis/libparmetis.a

# Define the required Fortran libraries, if you use C compiler to link
FLIBS	 	= 

# Define all the libraries
LIBS            = $(DSUPERLULIB) $(BLASLIB) $(PARMETISLIB) $(METISLIB)

#
#  The archiver and the flag(s) to use when building archive (library)
#  If your system has no ranlib, set RANLIB = echo.
#
ARCH         	= ar
ARCHFLAGS    	= cr
RANLIB       	= ranlib

############################################################################
# C compiler setup
CC           	= mpicc
# CFLAGS should be set to be the C flags that include optimization
CFLAGS          = -pipe -O2 -DPRNTlevel=2 -DDEBUGlevel=2
#
# NOOPTS should be set to be the C flags that turn off any optimization
NOOPTS		= 
############################################################################
# FORTRAN compiler setup
FORTRAN         = mpif77
F90FLAGS	= 
############################################################################
LOADER	        = mpif77
LOADOPTS	= 
############################################################################
#  C preprocessor defs for compilation (-DNoChange, -DAdd_, or -DUpCase)
#
#  Need follow the convention of how C calls a Fortran routine.
#
CDEFS        = -DAdd__ 
