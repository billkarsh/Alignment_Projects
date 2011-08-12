

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

# =====
# paths
# =====

PATH_JSON := /groups/scheffer/home/schefferl/json/jsoncpp-src-0.5.0/include

# =======
# sources
# =======

SOURCES_MYPROJ := \
	mos.cpp json_reader.cpp json_writer.cpp json_value.cpp

_SRC_UTL_ := \
	PipeFiles.cpp numerical_recipes.cpp mrc.cpp LinEqu.cpp

SOURCES_UTL := \
	$(PATH_UTL)/{$(subst $(space),$(comma),$(_SRC_UTL_))}

_SRC_GEN_ := \
	Cmdline.cpp File.cpp Disk.cpp ImageIO.cpp Maths.cpp Correlation.cpp \
	CPoint.cpp Geometry.cpp CTForm.cpp \
	Draw.cpp Memory.cpp

SOURCES_GEN :=  \
	$(PATH_GEN)/{$(subst $(space),$(comma),$(_SRC_GEN_))}

# =====
# build
# =====

mos:
	$(CC) $(CFLAGS) $(SYMBOLS_MYPROJ) \
$(OUTPUT) \
    \
$(INCLUDES_STD) \
    \
-I $(PATH_JSON) \
    \
$(SOURCES_MYPROJ) \
    \
$(SOURCES_UTL) \
    \
$(SOURCES_GEN) \
    \
$(OBJFILES_TIF) \
    \
$(STATICLIBS_FFT) \
    \
$(OTHERLIBS_STD) \
#### end $(CC) $@


