

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
	makemontages.cpp

_SRC_GEN_ := \
	Cmdline.cpp Disk.cpp File.cpp PipeFiles.cpp \
	CTileSet.cpp Maths.cpp Geometry.cpp CPoint.cpp \
	TAffine.cpp

SOURCES_GEN :=  \
	$(PATH_GEN)/{$(subst $(space),$(comma),$(_SRC_GEN_))}

# =====
# build
# =====

makemontages:
	$(CC) $(CFLAGS) $(SYMBOLS_MYPROJ) \
$(OUTPUT) \
    \
$(INCLUDES_STD) \
    \
$(SOURCES_MYPROJ) \
    \
$(SOURCES_GEN) \
    \
$(LINKS_STD) \
#### end $(CC) $@


