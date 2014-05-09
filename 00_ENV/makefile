
include $(ALN_LOCAL_MAKE_PATH)/aln_makefile_std_defs

1dirs =\
 1_CpyFMs\
 1_Cross_CarveBlocks\
 1_Cross_LowRes\
 1_Cross_Scaffold\
 1_Cross_ThisBlock\
 1_Cross_TopScripts\
 1_DelTiles\
 1_DMesh\
 1_EView\
 1_LSQi\
 1_LSQw\
 1_LsqErr\
 1_MakeIDB\
 1_MakeMontages\
 1_Mosaic_current\
 1_MRCSD1Lyr\
 1_Ptestx\
 1_Reformat\
 1_RemoveRefTiles\
 1_Scapeops\
 1_ShtMaker\
 1_Thumbs\
 1_Tiny\
 1_TopScripts\
 1_ViewSD\
 1_XView

2dirs =\
 2_EditXML\
 2_FFTomos\
 2_GraRan1Lyr\
 2_GraRanXML\
 2_GrayRanger\
 2_HEQ1Lyr\
 2_HEQ2Normal\
 2_HEQLayers\
 2_HEQVFY\
 2_HEQXML\
 2_Hist1\
 2_HistAll\
 2_LocalProj\
 2_ReverseTiles\
 2_RGBM1Lyr\
 2_RGBMerge\
 2_RGBMXML\
 2_SetApigChan\
 2_TFormTableEx\
 2_UniqueTiles\
 2_XMLAppend\
 2_XMLCopyOids\
 2_XMLCopyTF\
 2_XMLExists\
 2_XMLExtract\
 2_XMLGetTF\
 2_XMLInsert\
 2_XMLInterleave\
 2_XMLKillZ\
 2_XMLLocalPaths\
 2_XMLMissingZ\
 2_XMLRenum\
 2_XMLRevertTF\
 2_XMLRotTiles\
 2_XMLUnite\
 2_XMLZScale

all : xdirs = ${1dirs} ${2dirs}
all : build_dirs

all1 : xdirs = ${1dirs}
all1 : build_dirs

all2 : xdirs = ${2dirs}
all2 : build_dirs

clean : xdirs = ${1dirs} ${2dirs}
clean : clean_dirs clean_libs

build_dirs : .CHECK_ALL_LIBS
	@ \
	cd ..; \
	for i in $(xdirs); \
	do \
		cd $$i && $(MAKE) && cd ..; \
	done

clean_dirs :
	@ \
	cd ..; \
	for i in $(xdirs); \
	do \
		cd $$i && $(MAKE) clean && cd ..; \
	done

clean_libs:
	@ \
	cd ..; \
	for i in 0_BK 0_GEN 0_UTL 0_XML; \
	do \
		cd $$i && $(MAKE) clean && cd ..; \
	done

.PHONY : build_dirs clean_dirs clean_libs

