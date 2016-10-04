
include $(ALN_LOCAL_MAKE_PATH)/aln_makefile_std_defs

appname = local

files =\
 local.cpp

objs = ${files:.cpp=.o}

all : $(appname)

clean :
    rm -f *.o

$(appname) : .CHECK_GENLIB .DEV_MKDIR ${objs}
    $(CC) $(LFLAGS) ${objs} $(LINKS_STD) $(DEV_OUTPUT)

