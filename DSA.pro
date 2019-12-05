TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += \
        bigint.cpp \
        main.cpp \
        mpn.cpp \
        sha256.cpp

HEADERS += \
    bigint.h \
    mpn.h \
    sha256.h
