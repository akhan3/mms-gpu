CXX         := 	g++
CC         	:=	gcc
LINKER     	:= 	g++ -fPIC

#INCPATH     := -I./ $(MATLAB_INCLUDES)
#LIBPATH		:= $(MATLAB_LIBRARIES)
LIBRARIES	:= -lfreeimage

ifeq ($(dbg),1)
	COMMONFLAGS += -g
else
	COMMONFLAGS += -O3
endif

ifeq ($(prof),1)
	COMMONFLAGS += -pg
endif

CXXFLAGS    := -Wall -W $(INCPATH) $(COMMONFLAGS)
OBJS        := 	Box.o \
				Queue.o \
				Cmpx.o \
				Vector3.o \
				helper_functions.o \
				vector_functions.o \
				ode_functions.o \
				fmm_calc.o \
				main.o

TARGET    	:= main

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS) $(LIBPATH) $(LIBRARIES)

Box.o 				: Box.cpp Box.hpp
Queue.o				: Queue.cpp Queue.hpp
Cmpx.o				: Cmpx.cpp Cmpx.hpp
Vector3.o			: Vector3.cpp Vector3.hpp
helper_functions.o 	: helper_functions.cpp helper_functions.hpp
vector_functions.o 	: vector_functions.cpp vector_functions.hpp
ode_functions.o 	: ode_functions.cpp ode_functions.hpp
fmm_calc.o			: fmm_calc.cpp
main.o				: main.cpp

clean:
	rm -f $(OBJS) gmon.out
