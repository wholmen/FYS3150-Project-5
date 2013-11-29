TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt
LIBS += -larmadillo -lblas -llapack

SOURCES += main.cpp \
    pdesolver.cpp \
    functions.cpp \
    lib.cpp

HEADERS += \
    pdesolver.h \
    functions.h \
    lib.h

