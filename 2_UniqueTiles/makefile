
include $(ALN_LOCAL_MAKE_PATH)/aln_makefile_std_defs

appname = UniqueTiles

files =\
 UniqueTiles.cpp

objs = ${files:.cpp=.o}

all : $(appname)

clean :
    rm -f *.o

$(appname) : .CHECK_GENLIB ${objs}
    $(CC) $(LFLAGS) ${objs} $(LINKS_STD) $(OUTPUT)

