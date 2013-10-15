
include $(ALN_LOCAL_MAKE_PATH)/aln_makefile_std_defs

# =======
# defines
# =======

# -DCREATEMESH_WRITE_DEBUG_IMAGES:
#	Output 'Bogus3.tif' if can't find boundary left edge
#	Output 'Bogus4.tif' if isolated pixel on boundary

CFLAGS += -DCREATEMESH_WRITE_DEBUG_IMAGES

appname = ptest

files =\
 dmeshdriver.cpp\
 CGBL_dmesh.cpp\
 CThmUtil.cpp\
 ApproximateMatch_NoCR.cpp\
 ApproximateMatch.cpp\
 dmesh.cpp\
 RegionToRegionMap.cpp\
 CreateMesh.cpp\
 ImproveMesh.cpp\
 InSectionOverlap.cpp

objs = ${files:.cpp=.o}

all : $(appname)

clean :
	rm -f *.o

$(appname) : .CHECK_GENLIB ${objs}
	$(CC) $(LFLAGS) ${objs} $(LINKS_STD) $(OUTPUT)

