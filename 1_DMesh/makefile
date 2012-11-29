

# ==========
# ==========
# My Project
# ==========
# ==========

include $(ALN_LOCAL_MAKE_PATH)/aln_makefile_std_defs

# =======
# defines
# =======

# -DCREATEMESH_WRITE_DEBUG_IMAGES:
#	Output 'Bogus3.tif' if can't find boundary left edge
#	Output 'Bogus4.tif' if isolated pixel on boundary

SYMBOLS_MYPROJ := -DCREATEMESH_WRITE_DEBUG_IMAGES

# =======
# sources
# =======

SOURCES_MYPROJ := \
	dmeshdriver.cpp CGBL_dmesh.cpp CThmUtil.cpp \
	ApproximateMatch_NoCR.cpp ApproximateMatch.cpp \
	dmesh.cpp RegionToRegionMap.cpp CreateMesh.cpp \
	ImproveMesh.cpp InSectionOverlap.cpp

_SRC_UTL_ := \
	numerical_recipes.cpp mrc.cpp

SOURCES_UTL := \
	$(PATH_UTL)/{$(subst $(space),$(comma),$(_SRC_UTL_))}

_SRC_GEN_ := \
	Cmdline.cpp Disk.cpp File.cpp PipeFiles.cpp \
	ImageIO.cpp FoldMask.cpp Maths.cpp Correlation.cpp \
	CPoint.cpp Geometry.cpp CTForm.cpp CPixPair.cpp \
	CThmScan.cpp Metrics.cpp Inspect.cpp Timer.cpp \
	Memory.cpp Debug.cpp

SOURCES_GEN :=  \
	$(PATH_GEN)/{$(subst $(space),$(comma),$(_SRC_GEN_))}

# =====
# build
# =====

ptest:
	$(CC) $(CFLAGS) $(SYMBOLS_MYPROJ) \
$(OUTPUT) \
    \
$(INCLUDES_STD) \
    \
$(SOURCES_MYPROJ) \
    \
$(SOURCES_UTL) \
    \
$(SOURCES_GEN) \
    \
$(LINKS_STD) \
#### end $(CC) $@


