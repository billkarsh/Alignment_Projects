
include $(ALN_LOCAL_MAKE_PATH)/aln_makefile_std_defs

# =======
# defines
# =======

# -DTINYSTAT: Really lives in Maths.cpp, so we need to revisit this
# when doing foldmasks for real.
#

CFLAGS += -g -DTINYSTAT

appname = importMRC

files =\
 importMRC.cpp

objs = ${files:.cpp=.o}

all : $(appname)

clean :
    rm -f *.o

$(appname) : .CHECK_GENLIB .DEV_MKDIR ${objs}
    $(CC) $(LFLAGS) ${objs} $(LINKS_STD) $(DEV_OUTPUT)

