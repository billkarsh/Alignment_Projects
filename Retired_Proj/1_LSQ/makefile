
include $(ALN_LOCAL_MAKE_PATH)/aln_makefile_std_defs

appname = lsq

files =\
 lsq.cpp\
 lsq_Types.cpp\
 lsq_ReadPts.cpp\
 lsq_CNX.cpp\
 lsq_SML.cpp\
 lsq_DIR.cpp\
 lsq_MDL.cpp\
 lsq_MTrans.cpp\
 lsq_MSimlr.cpp\
 lsq_MAffine.cpp\
 lsq_MHmgphy.cpp

objs = ${files:.cpp=.o}

all : $(appname)

clean :
	rm -f *.o

$(appname) : .CHECK_GENLIB ${objs}
	$(CC) $(LFLAGS) ${objs} $(LINKS_STD) $(OUTPUT)

