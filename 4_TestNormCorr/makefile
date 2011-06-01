

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

SYMBOLS_MYPROJ :=

# =======
# sources
# =======

SOURCES_MYPROJ := \
	testnormcorr.cpp

_SRC_GEN_ := \
	Maths.cpp Correlation.cpp \
	Geometry.cpp Timer.cpp

SOURCES_GEN :=  \
	$(PATH_GEN)/{$(subst $(space),$(comma),$(_SRC_GEN_))}

# =====
# build
# =====

testnormcorr:
	$(shell mkdir -p $(PATH_OUT)) \
	$(CC) $(CFLAGS) $(SYMBOLS_MYPROJ) \
$(OUTPUT) \
    \
$(INCLUDES_STD) \
    \
$(SOURCES_MYPROJ) \
    \
$(SOURCES_GEN) \
    \
$(OBJFILES_TIF) \
    \
$(STATICLIBS_FFT) \
    \
$(OTHERLIBS_STD) \
#### end $(CC) $@


