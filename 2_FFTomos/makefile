
include $(ALN_LOCAL_MAKE_PATH)/aln_makefile_std_defs

CFLAGS += -I$(PATH_BK)

appname = fftomos

files =\
 FFTomos.cpp\
 CMask.cpp

objs = ${files:.cpp=.o}

bk_files =\
 BK_BMAP.c\
 BK_BMAP_CONVERT.c\
 BK_BMAP_SET.c\
 BK_BMP.c\
 BK_GEOM.c\
 BK_HST.c\
 BK_MEM.c\
 BK_RGN.c\
 BK_SORT.c\
 BK_SORT_FP32.c\
 BK_STAT.c\
 BK_SUM.c

bk_src = $(PATH_BK)/{$(subst $(space),$(comma),$(bk_files))}

all : $(appname)

clean :
    rm -f *.o

$(appname) : .CHECK_GENLIB ${objs}
    $(CC) $(CFLAGS) $(bk_src) $(LFLAGS) ${objs} $(LINKS_STD) $(OUTPUT)

