

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
	FFTomos.cpp CMask.cpp

_SRC_UTL_ := \
	numerical_recipes.cpp mrc.cpp

SOURCES_UTL := \
	$(PATH_UTL)/{$(subst $(space),$(comma),$(_SRC_UTL_))}

_SRC_GEN_ := \
	Cmdline.cpp Disk.cpp File.cpp \
	ImageIO.cpp Maths.cpp TAffine.cpp

SOURCES_GEN :=  \
	$(PATH_GEN)/{$(subst $(space),$(comma),$(_SRC_GEN_))}

_SRC_BK_ := \
	BK_BMAP.c BK_BMAP_CONVERT.c BK_BMAP_SET.c \
	BK_BMP.c BK_GEOM.c BK_HST.c \
	BK_MEM.c BK_RGN.c BK_SORT.c \
	BK_SORT_FP32.c BK_STAT.c BK_SUM.c

SOURCES_BK :=  \
	$(PATH_BK)/{$(subst $(space),$(comma),$(_SRC_BK_))}

# =====
# build
# =====

fftomos:
	$(CC) $(CFLAGS) $(SYMBOLS_MYPROJ) \
$(OUTPUT) \
    \
$(INCLUDES_STD) -I$(PATH_BK) \
    \
$(SOURCES_MYPROJ) \
    \
$(SOURCES_UTL) \
    \
$(SOURCES_GEN) \
    \
$(SOURCES_BK) \
    \
$(LINKS_STD) \
#### end $(CC) $@


