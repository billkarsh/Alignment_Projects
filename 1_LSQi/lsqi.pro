DEPENDPATH += $$PWD
INCLUDEPATH += $$PWD

INCLUDEPATH += ../0_GEN
include(../0_GEN/gen.pro)

HEADERS += \
    $$PWD/lsq_Layers.h

SOURCES += \
    $$PWD/lsq.cpp \
    $$PWD/lsq_Layers.cpp

