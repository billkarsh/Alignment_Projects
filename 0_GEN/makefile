
include $(ALN_LOCAL_MAKE_PATH)/aln_makefile_std_defs

libname = genlib

files =\
 CAffineLens.cpp\
 CCorrImages.cpp\
 CCropMask.cpp\
 Cmdline.cpp\
 Correlation.cpp\
 CPicBase.cpp\
 CPixPair.cpp\
 CPoint.cpp\
 CRegexID.cpp\
 CRigid.cpp\
 CTemplate.cpp\
 CThmScan.cpp\
 CTileSet.cpp\
 CTileSet_Scape.cpp\
 Debug.cpp\
 Disk.cpp\
 Draw.cpp\
 EZThreads.cpp\
 File.cpp\
 FoldMask.cpp\
 Geometry.cpp\
 ImageIO.cpp\
 Inspect.cpp\
 LinEqu.cpp\
 Maths.cpp\
 Memory.cpp\
 Metrics.cpp\
 PipeFiles.cpp\
 Scape.cpp\
 TAffine.cpp\
 Tform_Array.cpp\
 THmgphy.cpp\
 Timer.cpp\
 TrakEM2_UTL.cpp

objs = ${files:.cpp=.o}

all : $(libname).a

clean :
	rm -f *.{a,o}

$(libname).a : ${objs}
	ar rvs $(libname).a ${objs}

