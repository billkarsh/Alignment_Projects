

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
	fft.cpp

# =====
# build
# =====

fft:
	$(shell mkdir -p $(PATH_OUT)) \
	$(CC) $(CFLAGS) $(SYMBOLS_MYPROJ) \
$(OUTPUT) \
    \
$(INCLUDES_STD) \
    \
$(SOURCES_MYPROJ) \
    \
$(STATICLIBS_FFT) \
    \
$(OTHERLIBS_STD) \
#### end $(CC) $@


