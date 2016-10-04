
include $(ALN_LOCAL_MAKE_PATH)/aln_makefile_std_defs

#PATH_JSON := /groups/scheffer/home/schefferl/json/jsoncpp-src-0.5.0/include
PATH_JSON := json-include
CFLAGS    += -I$(PATH_JSON)

appname = mos

files =\
 mos.cpp\
 json_reader.cpp\
 json_writer.cpp\
 json_value.cpp

objs = ${files:.cpp=.o}

all : $(appname)

clean :
    rm -f *.o

$(appname) : .CHECK_GENLIB ${objs}
    $(CC) $(LFLAGS) ${objs} $(LINKS_STD) $(OUTPUT)

