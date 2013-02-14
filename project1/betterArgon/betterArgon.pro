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

LIBS += -larmadillo -llapack -lblas

QMAKE_CXXFLAGS -= -O
QMAKE_CXXFLAGS -= -O2

QMAKE_CXXFLAGS += -O3
