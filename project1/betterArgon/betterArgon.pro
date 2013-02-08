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
