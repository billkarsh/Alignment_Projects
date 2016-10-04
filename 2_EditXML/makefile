
include $(ALN_LOCAL_MAKE_PATH)/aln_makefile_std_defs

appname = EditXML

files =\
 EditXML.cpp

objs = ${files:.cpp=.o}

all : $(appname)

clean :
	rm -f *.o

$(appname) : .CHECK_GENLIB ${objs}
	$(CC) $(LFLAGS) ${objs} $(LINKS_STD) $(OUTPUT)

