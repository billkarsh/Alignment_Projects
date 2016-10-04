
include $(ALN_LOCAL_MAKE_PATH)/aln_makefile_std_defs

appname = cross_topscripts

files =\
 cross_topscripts.cpp

objs = ${files:.cpp=.o}

all : $(appname)

clean :
    rm -f *.o

$(appname) : .CHECK_GENLIB ${objs}
    $(CC) $(LFLAGS) ${objs} $(LINKS_STD) $(OUTPUT)

