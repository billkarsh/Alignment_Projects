
include $(ALN_LOCAL_MAKE_PATH)/aln_makefile_std_defs

appname = thumbs

files =\
 ThumbsDriver.cpp\
 CGBL_Thumbs.cpp\
 ../1_DMesh/CThmUtil.cpp\
 Thumbs_NoCR.cpp\
 Thumbs.cpp

objs = ${files:.cpp=.o}

all : $(appname)

clean :
	rm -f *.o

$(appname) : .CHECK_GENLIB ${objs}
	$(CC) $(LFLAGS) ${objs} $(LINKS_STD) $(OUTPUT)

