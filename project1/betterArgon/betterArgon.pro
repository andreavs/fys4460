TEMPLATE = app
CONFIG += console
CONFIG -= qt

SOURCES += main.cpp \
    atom.cpp \
    system.cpp \
    cell.cpp

HEADERS += \
    atom.h \
    system.h \
    cell.h

LIBS += -larmadillo -llapack -lblas -fopenmp

COMMON_CXXFLAGS = -std=c++0x -fopenmp
QMAKE_CXXFLAGS += $$COMMON_CXXFLAGS
QMAKE_CXXFLAGS_RELEASE += $$COMMON_CXXFLAGS
QMAKE_CXXFLAGS_DEBUG += $$COMMON_CXXFLAGS

release {
DEFINES += ARMA_NO_DEBUG
QMAKE_LFLAGS -= -O1
QMAKE_LFLAGS += -O3
QMAKE_LFLAGS_RELEASE -= -O1
QMAKE_LFLAGS_RELEASE += -O3
QMAKE_CXXFLAGS -= -O2
QMAKE_CXXFLAGS += -O3
QMAKE_CXXFLAGS_RELEASE -= -O2
QMAKE_CXXFLAGS_RELEASE += -O3
}
