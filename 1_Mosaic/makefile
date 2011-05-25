

# ==========
# ==========
# My Project
# ==========
# ==========

include ../00_ENV/makefile_std_defs

# =======
# defines
# =======

SYMBOLS_MYPROJ :=

# =======
# sources
# =======

SOURCES_MYPROJ := \
	mos.cpp

_SRC_UTL_ := \
	numerical_recipes.cpp mrc.cpp LinEqu.cpp

SOURCES_UTL := \
	$(PATH_UTL)/{$(subst $(space),$(comma),$(_SRC_UTL_))}

_SRC_GEN_ := \
	File.cpp \
	ImageIO.cpp Maths.cpp Correlation.cpp \
	CPoint.cpp Geometry.cpp CTForm.cpp \
	Draw.cpp

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


