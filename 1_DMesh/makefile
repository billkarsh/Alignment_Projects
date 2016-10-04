
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
 ApproximateMatch_NoCR.cpp\
 ApproximateMatch.cpp\
 CGBL_dmesh.cpp\
 CreateMesh.cpp\
 CThmUtil.cpp\
 dmesh.cpp\
 dmeshdriver.cpp\
 ImproveMesh.cpp\
 InSectionOverlap.cpp\
 janelia.cpp\
 RegionToRegionMap.cpp

objs = ${files:.cpp=.o}

all : $(appname)

clean :
	rm -f *.o

$(appname) : .CHECK_GENLIB ${objs}
	$(CC) $(LFLAGS) ${objs} $(LINKS_STD) $(OUTPUT)

