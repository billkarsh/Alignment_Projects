

# ==========
# ==========
# My Project
# ==========
# ==========

include $(ALN_LOCAL_MAKE_PATH)/aln_makefile_std_defs

# =======
# defines
# =======

SYMBOLS_MYPROJ := -DTIXML_USE_STL

# =======
# sources
# =======

SOURCES_MYPROJ := \
	cross_carveblocks.cpp

_SRC_XML_ := \
	tinyxml.cpp tinyxmlerror.cpp tinyxmlparser.cpp

SOURCES_XML := \
	$(PATH_XML)/{$(subst $(space),$(comma),$(_SRC_XML_))}

_SRC_GEN_ := \
	Cmdline.cpp Disk.cpp File.cpp PipeFiles.cpp \
	TrakEM2_UTL.cpp CTileSet.cpp Geometry.cpp TAffine.cpp

SOURCES_GEN :=  \
	$(PATH_GEN)/{$(subst $(space),$(comma),$(_SRC_GEN_))}

# =====
# build
# =====

cross_carveblocks:
	$(CC) $(CFLAGS) $(SYMBOLS_MYPROJ) \
$(OUTPUT) \
    \
$(INCLUDES_STD) \
    \
$(SOURCES_MYPROJ) \
    \
$(SOURCES_XML) \
    \
$(SOURCES_GEN) \
    \
$(LINKS_STD) \
#### end $(CC) $@


