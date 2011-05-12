

# ==========
# ==========
# My Project
# ==========
# ==========

include ../00_ENV/makefile_std_defs

# =======
# defines
# =======

SYMBOLS_MYPROJ := -DTIXML_USE_STL

# =======
# sources
# =======

SOURCES_MYPROJ := \
	XMLMissingZ.cpp

_SRC_XML_ := \
	tinyxml.cpp tinyxmlerror.cpp tinyxmlparser.cpp

SOURCES_XML := \
	$(PATH_XML)/{$(subst $(space),$(comma),$(_SRC_XML_))}

_SRC_GEN_ := \
	Cmdline.cpp File.cpp

SOURCES_GEN :=  \
	$(PATH_GEN)/{$(subst $(space),$(comma),$(_SRC_GEN_))}

# =====
# build
# =====

XMLMissingZ:
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
$(OBJFILES_TIF) \
    \
$(STATICLIBS_FFT) \
    \
$(OTHERLIBS_STD) \
#### end $(CC) $@


