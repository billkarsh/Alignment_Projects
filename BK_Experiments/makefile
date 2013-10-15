
include $(ALN_LOCAL_MAKE_PATH)/aln_makefile_std_defs

DEGUG = -g

targets =\
 davisubset\
 diff\
 fixcoords\
 junk\
 linesolap\
 pngtest\
 temcoordfix\
 test\
 vectors

all : $(targets)

davisubset : davisubset.o .CHECK_GENLIB
	$(CC) $(CFLAGS) $< $(LFLAGS) $(LINKS_STD) $(OUTPUT)

diff : diff.o .CHECK_GENLIB
	$(CC) $(CFLAGS) $< $(LFLAGS) $(LINKS_STD) $(OUTPUT)

fixcoords : fixcoords.o .CHECK_GENLIB
	$(CC) $(CFLAGS) $< $(LFLAGS) $(LINKS_STD) $(OUTPUT)

junk : junk.o .CHECK_GENLIB
	$(CC) $(CFLAGS) $< $(LFLAGS) $(LINKS_STD) $(OUTPUT)

linesolap : linesolap.o .CHECK_GENLIB
	$(CC) $(CFLAGS) $< $(LFLAGS) $(LINKS_STD) $(OUTPUT)

pngtest : pngtest.o .CHECK_GENLIB
	$(CC) $(CFLAGS) $(DEBUG) $< $(LFLAGS) $(LINKS_STD) $(OUTPUT)

temcoordfix : temcoordfix.o .CHECK_GENLIB
	$(CC) $(CFLAGS) $< $(LFLAGS) $(LINKS_STD) $(OUTPUT)

test : test.o .CHECK_GENLIB
	$(CC) $(CFLAGS) $< $(LFLAGS) $(LINKS_STD) $(OUTPUT)

vectors : vectors.o .CHECK_GENLIB
	$(CC) $(CFLAGS) $(DEBUG) $< $(LFLAGS) $(LINKS_STD) $(OUTPUT)

clean :
	rm -f *.o

.PHONY : $(targets)

