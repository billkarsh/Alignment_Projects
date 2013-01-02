

# ==========
# ==========
# My Project
# ==========
# ==========

include $(ALN_LOCAL_MAKE_PATH)/aln_makefile_std_defs

# =======
# defines
# =======

SYMBOLS_MYPROJ :=

# =======
# sources
# =======

SOURCES_MYPROJ := \
	ThumbsDriver.cpp CGBL_Thumbs.cpp ../1_DMesh/CThmUtil.cpp \
	Thumbs_NoCR.cpp Thumbs.cpp

_SRC_UTL_ := \
	numerical_recipes.cpp mrc.cpp

SOURCES_UTL := \
	$(PATH_UTL)/{$(subst $(space),$(comma),$(_SRC_UTL_))}

_SRC_GEN_ := \
	Cmdline.cpp Disk.cpp File.cpp PipeFiles.cpp \
	CLens.cpp ImageIO.cpp FoldMask.cpp Maths.cpp \
	Correlation.cpp CPoint.cpp Geometry.cpp CTForm.cpp \
	CPixPair.cpp CThmScan.cpp Timer.cpp Memory.cpp \
	Debug.cpp

SOURCES_GEN :=  \
	$(PATH_GEN)/{$(subst $(space),$(comma),$(_SRC_GEN_))}

# =====
# build
# =====

thumbs:
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


