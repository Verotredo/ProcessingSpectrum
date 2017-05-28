#-------------------------------------------------
#
# Project created by QtCreator 2017-05-18T20:45:02
#
#-------------------------------------------------

QT       += core gui charts

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = ProcessingPlots
TEMPLATE = app

# The following define makes your compiler emit warnings if you use
# any feature of Qt which as been marked as deprecated (the exact warnings
# depend on your compiler). Please consult the documentation of the
# deprecated API in order to know how to port your code away from it.
DEFINES += QT_DEPRECATED_WARNINGS

# You can also make your code fail to compile if you use deprecated APIs.
# In order to do so, uncomment the following line.
# You can also select to disable deprecated APIs only up to a certain version of Qt.
#DEFINES += QT_DISABLE_DEPRECATED_BEFORE=0x060000    # disables all the APIs deprecated before Qt 6.0.0


SOURCES += main.cpp\
    processplot.cpp \
    matvector.cpp \
    plotwidget.cpp \
    cpp/src/alglibinternal.cpp \
    cpp/src/alglibmisc.cpp \
    cpp/src/ap.cpp \
    cpp/src/dataanalysis.cpp \
    cpp/src/diffequations.cpp \
    cpp/src/fasttransforms.cpp \
    cpp/src/integration.cpp \
    cpp/src/interpolation.cpp \
    cpp/src/linalg.cpp \
    cpp/src/optimization.cpp \
    cpp/src/solvers.cpp \
    cpp/src/specialfunctions.cpp \
    cpp/src/statistics.cpp

HEADERS  += \
    processplot.h \
    matvector.h \
    plotwidget.h \
    cpp/src/alglibinternal.h \
    cpp/src/alglibmisc.h \
    cpp/src/ap.h \
    cpp/src/dataanalysis.h \
    cpp/src/diffequations.h \
    cpp/src/fasttransforms.h \
    cpp/src/integration.h \
    cpp/src/interpolation.h \
    cpp/src/linalg.h \
    cpp/src/optimization.h \
    cpp/src/solvers.h \
    cpp/src/specialfunctions.h \
    cpp/src/statistics.h \
    cpp/src/stdafx.h

FORMS    +=