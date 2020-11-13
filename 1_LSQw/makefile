
include $(ALN_LOCAL_MAKE_PATH)/aln_makefile_std_defs

ifeq ($(USING_MPI), 1)
# If Valgrind (run as valgrind --leak-check=yes <mycmdline>)
#CFLAGS = -g -O0
#LFLAGS =

# Compile mpi-based apps with CC := mpicxx (gnu), or mpiicpc (intel)
# and add -mt_mpi switch to LFLAGS
CC := mpicxx
LFLAGS += -mt_mpi
endif

appname = lsqw

files =\
 lsqw.cpp\
 ../1_LSQi/lsq_Layers.cpp\
 lsq_Bounds.cpp\
 lsq_Dropout.cpp\
 lsq_Error.cpp\
 lsq_Globals.cpp\
 lsq_LoadPoints.cpp\
 lsq_Magnitude.cpp\
 lsq_MPI.cpp\
 lsq_Solve.cpp\
 lsq_Split.cpp\
 lsq_Untwist.cpp\
 lsq_XArray.cpp

objs = ${files:.cpp=.o}

all : $(appname)

clean :
    rm -f *.o

$(appname) : .CHECK_GENLIB ${objs}
    $(CC) $(LFLAGS) ${objs} $(LINKS_STD) $(OUTPUT)

