CXX         := 	g++
CC         	:=	gcc
LINKER     	:= 	g++ -fPIC

#INCPATH     := -I./ $(MATLAB_INCLUDES)
#LIBPATH		:= $(MATLAB_LIBRARIES)
LIBRARIES	:= -lfreeimage
ifeq ($(dbg),1)
	#COMMONFLAGS += -g -pg
	COMMONFLAGS += -g
else
	#COMMONFLAGS += -O3 -pg
	COMMONFLAGS += -O3
endif
CXXFLAGS    := -Wall -W $(INCPATH) $(COMMONFLAGS)
OBJS        := main.o Box.o Queue.o Cmpx.o fmm_calc.o helper_functions.o
TARGET    	:= main

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS) $(LIBPATH) $(LIBRARIES)

main.o		: main.cpp Box.hpp Queue.hpp helper_functions.hpp
Box.o 		: Box.cpp Box.hpp
Queue.o		: Queue.cpp Queue.hpp
Cmpx.o		: Cmpx.cpp Cmpx.hpp
fmm_calc.o	: fmm_calc.cpp Box.hpp Queue.hpp Cmpx.hpp helper_functions.hpp
helper_functions.o : helper_functions.cpp helper_functions.hpp

clean:
	rm -f $(OBJS) $(TARGET) gmon.out
