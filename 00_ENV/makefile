
1dirs =\
 1_Align\
 1_Clean_ThmPairs\
 1_CpyFMs\
 1_Cross_CarveBlocks\
 1_Cross_LowRes\
 1_Cross_LowToHiRes\
 1_Cross_ThisBlock\
 1_Cross_TopScripts\
 1_DelTiles\
 1_DMesh\
 1_LSQ\
 1_Lsqerr\
 1_MakeALN\
 1_MakeIDB\
 1_MakeMontages\
 1_Make_ldir\
 1_Mosaic_current\
 1_MRCSD1Lyr\
 1_Ptestx\
 1_Reformat\
 1_RemoveRefTiles\
 1_Scapeops\
 1_Script\
 1_ShtMaker\
 1_Thumbs\
 1_Tiny\
 1_TopScripts\
 1_UpdateMontages\
 1_ViewSD

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
 2_XMLLocalPaths\
 2_XMLMissingZ\
 2_XMLRenum\
 2_XMLRevertTF\
 2_XMLRotTiles\
 2_XMLUnite\
 2_XMLZScale

4dirs =\
 4_Cmap\
 4_Connect\
 4_FFT\
 4_Fold\
 4_Import\
 4_Local\
 4_Master\
 4_Qa\
 4_Qafm\
 4_Qcheck\
 4_Sla\
 4_TestNormCorr\
 4_Time_FFT\
 4_Two\
 4_Two32\
 4_Unfold\
 4_Unmap

all: xdirs = ${1dirs} ${2dirs}
all: build_loop

all1: xdirs = ${1dirs}
all1: build_loop

all2: xdirs = ${2dirs}
all2: build_loop

all4: xdirs = ${4dirs}
all4: build_loop

build_loop:
	@ \
	cd ..; \
	for i in $(xdirs); \
	do \
		cd $$i && $(MAKE) && cd ..; \
	done



