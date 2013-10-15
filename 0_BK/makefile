
include $(ALN_LOCAL_MAKE_PATH)/aln_makefile_std_defs

libname =

files =

objs = ${files:.cpp=.o}

all : 

clean :
	rm -f *.{a,o}

$(libname).a : ${objs}
	ar rvs $(libname).a ${objs}

