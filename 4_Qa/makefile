

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

SYMBOLS_MYPROJ := -DTINYSTAT

# =======
# sources
# =======

SOURCES_MYPROJ := \
	qa.cpp

_SRC_UTL_ := \
	numerical_recipes.cpp mrc.cpp

SOURCES_UTL := \
	$(PATH_UTL)/{$(subst $(space),$(comma),$(_SRC_UTL_))}

_SRC_GEN_ := \
	File.cpp ImageIO.cpp Maths.cpp \
	Correlation.cpp Geometry.cpp Debug.cpp

SOURCES_GEN :=  \
	$(PATH_GEN)/{$(subst $(space),$(comma),$(_SRC_GEN_))}

# =====
# build
# =====

qa:
	$(shell mkdir -p $(PATH_OUT)) \
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


