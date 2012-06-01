

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
	updatemontages.cpp

_SRC_XML_ := \
	tinyxml.cpp tinyxmlerror.cpp tinyxmlparser.cpp

SOURCES_XML := \
	$(PATH_XML)/{$(subst $(space),$(comma),$(_SRC_XML_))}

_SRC_UTL_ := \
	PipeFiles.cpp

SOURCES_UTL := \
	$(PATH_UTL)/$(_SRC_UTL_)

_SRC_GEN_ := \
	Cmdline.cpp CRegexID.cpp Disk.cpp File.cpp \
	TrakEM2_UTL.cpp CTForm.cpp

SOURCES_GEN :=  \
	$(PATH_GEN)/{$(subst $(space),$(comma),$(_SRC_GEN_))}

# =====
# build
# =====

updatemontages:
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


