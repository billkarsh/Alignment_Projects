DEPENDPATH += $$PWD
INCLUDEPATH += $$PWD

INCLUDEPATH += ../0_GEN
include(../0_GEN/gen.pro)

INCLUDEPATH += ../1_LSQi
include(../1_LSQi/lsqi.pro)

HEADERS += \
    $$PWD/lsq_Bounds.h \
    $$PWD/lsq_Dropout.h \
    $$PWD/lsq_Error.h \
    $$PWD/lsq_Globals.h \
    $$PWD/lsq_LoadPoints.h \
    $$PWD/lsq_Magnitude.h \
    $$PWD/lsq_MPI.h \
    $$PWD/lsq_Solve.h \
    $$PWD/lsq_Split.h \
    $$PWD/lsq_Untwist.h \
    $$PWD/lsq_XArray.h

SOURCES += \
    $$PWD/lsq_Bounds.cpp \
    $$PWD/lsq_Dropout.cpp \
    $$PWD/lsq_Error.cpp \
    $$PWD/lsq_Globals.cpp \
    $$PWD/lsq_LoadPoints.cpp \
    $$PWD/lsq_Magnitude.cpp \
    $$PWD/lsq_MPI.cpp \
    $$PWD/lsq_Solve.cpp \
    $$PWD/lsq_Split.cpp \
    $$PWD/lsq_Untwist.cpp \
    $$PWD/lsqw.cpp \
    $$PWD/lsq_XArray.cpp

