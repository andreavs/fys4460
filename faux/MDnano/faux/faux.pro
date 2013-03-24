TEMPLATE = app
CONFIG += console
CONFIG -= qt

SOURCES += main.cpp \
    system.cpp \
    cell.cpp \
    atom.cpp \
    simulation.cpp \
    statisticscalculator.cpp \
    configreader.cpp

HEADERS += \
    system.h \
    cell.h \
    atom.h \
    simulation.h \
    statisticscalculator.h \
    configreader.h

LIBS += -larmadillo -llapack -lblas -fopenmp -lconfig++

#LIBS += /opt/intel/lib/intel64/libomp5.so

#INCLUDEPATH += /opt/intel/lib/intel64

COMMON_CXXFLAGS = -std=c++0x -fopenmp
QMAKE_CXXFLAGS += $$COMMON_CXXFLAGS
QMAKE_CXXFLAGS_RELEASE += $$COMMON_CXXFLAGS
QMAKE_CXXFLAGS_DEBUG += $$COMMON_CXXFLAGS
QMAKE_LFLAGS += -openmp

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
