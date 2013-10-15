
include $(ALN_LOCAL_MAKE_PATH)/aln_makefile_std_defs

libname = utllib

files =\
 jh6.cpp\
 ls_svd.cpp\
 numerical_recipes.cpp\
 mrc.cpp

objs = ${files:.cpp=.o}

all : $(libname).a

clean :
	rm -f *.{a,o}

$(libname).a : ${objs}
	ar rvs $(libname).a ${objs}

