

# ==========
# ==========
# My Project
# ==========
# ==========

include $(ALN_LOCAL_MAKE_PATH)/aln_makefile_std_defs

# =====
# paths
# =====

PATH_OUT := $(PATH_OUT)/compiletest

# =======
# defines
# =======

SYMBOLS_MYPROJ := -DTIXML_USE_STL

# =======
# sources
# =======

SOURCES_MYPROJ := \
	master.cpp

_SRC_XML_ := \
	tinyxml.cpp tinyxmlerror.cpp tinyxmlparser.cpp

SOURCES_XML := \
	$(PATH_XML)/{$(subst $(space),$(comma),$(_SRC_XML_))}

_SRC_UTL_ := \
	jh6.cpp ls_svd.cpp \
	numerical_recipes.cpp mrc.cpp

SOURCES_UTL := \
	$(PATH_UTL)/{$(subst $(space),$(comma),$(_SRC_UTL_))}

_SRC_GEN_ := \
	File.cpp ImageIO.cpp Maths.cpp Correlation.cpp \
	CPoint.cpp Geometry.cpp TAffine.cpp CPicBase.cpp \
	CTemplate.cpp CCorrImages.cpp Debug.cpp

SOURCES_GEN :=  \
	$(PATH_GEN)/{$(subst $(space),$(comma),$(_SRC_GEN_))}

# =====
# build
# =====

master:
	$(shell mkdir -p $(PATH_OUT)) \
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


