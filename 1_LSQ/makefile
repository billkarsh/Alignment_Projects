

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
	lsq.cpp lsq_Types.cpp lsq_ReadPts.cpp \
	lsq_CNX.cpp lsq_SML.cpp lsq_DIR.cpp \
	lsq_MDL.cpp lsq_MTrans.cpp lsq_MSimlr.cpp \
	lsq_MAffine.cpp lsq_MHmgphy.cpp

_SRC_XML_ := \
	tinyxml.cpp tinyxmlerror.cpp tinyxmlparser.cpp

SOURCES_XML := \
	$(PATH_XML)/{$(subst $(space),$(comma),$(_SRC_XML_))}

_SRC_UTL_ := \
	numerical_recipes.cpp mrc.cpp

SOURCES_UTL := \
	$(PATH_UTL)/{$(subst $(space),$(comma),$(_SRC_UTL_))}

_SRC_GEN_ := \
	Cmdline.cpp CRegexID.cpp Disk.cpp File.cpp \
	PipeFiles.cpp CAffineLens.cpp TrakEM2_UTL.cpp \
	ImageIO.cpp Maths.cpp CPoint.cpp TAffine.cpp \
	THmgphy.cpp LinEqu.cpp Memory.cpp

SOURCES_GEN :=  \
	$(PATH_GEN)/{$(subst $(space),$(comma),$(_SRC_GEN_))}

# =====
# build
# =====

lsq:
	$(CC) $(CFLAGS) $(SYMBOLS_MYPROJ) \
$(OUTPUT) \
    \
$(INCLUDES_STD) \
    \
$(SOURCES_MYPROJ) \
    \
$(SOURCES_XML) \
    \
$(SOURCES_UTL) \
    \
$(SOURCES_GEN) \
    \
$(LINKS_STD) \
#### end $(CC) $@


