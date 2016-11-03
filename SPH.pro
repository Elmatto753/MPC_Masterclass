TARGET=SPHFluid
CONFIG+=c++11
CONFIG-= x86_64 app_bundle

OBJECTS_DIR=obj
MOC_DIR=moc
QT+= opengl core gui

INCLUDEPATH+=./include \
            /usr/local/include

SOURCES+=$$PWD/src/sph_main.cpp \
         $$PWD/src/sph_system.cpp \
         $$PWD/src/sph_timer.cpp \
         $$PWD/src/sph_phase.cpp \
         $$PWD/src/sph_particle.cpp

HEADERS+=$$PWD/include/sph_data.h \
         $$PWD/include/sph_header.h \
         $$PWD/include/sph_system.h \
         $$PWD/include/sph_timer.h \
         $$PWD/include/sph_type.h \
         $$PWD/include/sph_phase.h \
         $$PWD/include/sph_particle.h

OTHER_FILES+=$$PWD/Shader/shader.fs \
             $$PWD/Shader/shader.vs

DESTDIR=./

macx:LIBS+= -framework OpenGL

linux: LIBS+= -lGLU -lGLEW -lglut
