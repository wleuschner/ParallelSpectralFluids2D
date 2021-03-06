#-------------------------------------------------
#
# Project created by QtCreator 2018-09-14T01:44:02
#
#-------------------------------------------------

QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = PSF2D
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
        mainwindow.cpp \
    dec.cpp \
    mesh2d.cpp \
    vertex2d.cpp \
    edge2d.cpp \
    face2d.cpp \
    voxelviswidget.cpp \
    flowviswidget.cpp \
    imageviswidget.cpp \
    abstractsolver.cpp \
    decmesh2d.cpp \
    spectralfluidssolver2d.cpp \
    spectralfluidssolver2domp.cpp \
    abstractflowviswidget.cpp \
    spectralfluidssolver2dcl.cpp

HEADERS  += mainwindow.h \
    dec.h \
    mesh2d.h \
    vertex2d.h \
    edge2d.h \
    face2d.h \
    voxelviswidget.h \
    flowviswidget.h \
    imageviswidget.h \
    abstractsolver.h \
    decmesh2d.h \
    spectralfluidssolver2d.h \
    spectralfluidssolver2domp.h \
    abstractflowviswidget.h \
    spectralfluidssolver2dcl.h \
    gridenums.h

LIBS += -lOpenCL -larpack -lGLU
FORMS    += mainwindow.ui
DEFINES = GLM_ENABLE_EXPERIMENTAL GLM_FORCE_SSE2 VIENNACL_HAVE_EIGEN

QMAKE_CXXFLAGS += -std=c++17

QMAKE_CXXFLAGS_RELEASE -= -O2
QMAKE_CXXFLAGS_RELEASE += -O3 -march=native

QMAKE_LFLAGS_RELEASE -= -O1
QMAKE_LFLAGS_RELEASE += -O3

QMAKE_CXXFLAGS+= -fopenmp
QMAKE_LFLAGS +=  -fopenmp

DISTFILES += \
    flow.fs \
    flow.vs
