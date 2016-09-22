DEPENDPATH += $$PWD
INCLUDEPATH += $$PWD

INCLUDEPATH += ../0_GEN
include(../0_GEN/gen.pro)

HEADERS += \
    $$PWD/ApproximateMatch.h \
    $$PWD/CGBL_dmesh.h \
    $$PWD/CreateMesh.h \
    $$PWD/CThmUtil.h \
    $$PWD/dmesh.h \
    $$PWD/ImproveMesh.h \
    $$PWD/InSectionOverlap.h \
    $$PWD/RegionToRegionMap.h

SOURCES += \
    $$PWD/ApproximateMatch.cpp \
    $$PWD/ApproximateMatch_NoCR.cpp \
    $$PWD/CGBL_dmesh.cpp \
    $$PWD/CreateMesh.cpp \
    $$PWD/CThmUtil.cpp \
    $$PWD/dmesh.cpp \
    $$PWD/dmeshdriver.cpp \
    $$PWD/dmesh_unused.cpp \
    $$PWD/ImproveMesh.cpp \
    $$PWD/InSectionOverlap.cpp \
    $$PWD/RegionToRegionMap.cpp

