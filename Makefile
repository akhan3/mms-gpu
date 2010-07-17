#CXX         := 	g++-4.3
#CC         	:=	gcc-4.3
#LINKER     	:= 	g++-4.3 -fPIC
CXX         := 	g++
CC         	:=	gcc
LINKER     	:= 	g++ -fPIC

# My includes and libraries
INCLUDES   	+= -I. $(MATLAB_INCLUDES)
LIBRARIES 	+= $(MATLAB_LIBRARIES)

ifeq ($(dbg),1)
	CXXFLAGS 	+= 	-g -D__DEBUG__
	OBJDIR		+= 	obj/debug/
	BINDIR		+= 	bin/debug/
else
	CXXFLAGS 	+= 	-O3
	OBJDIR		+= 	obj/release/
	BINDIR		+= 	bin/release/
endif

ifeq ($(verbose),1)
	VERBOSE	:=
else
	VERBOSE	:=	@
endif

COMMONFLAGS		+=	-Wall -W $(INCLUDES)
CXXFLAGS    	+= 	$(COMMONFLAGS)

OBJS	:=	$(OBJDIR)quadtree.o
				
TARGET	:= 	$(BINDIR)quadtree

# ==================================================================
# Rules, target and dependencies
# ==================================================================

$(TARGET):	compile create_bin_dir
	$(VERBOSE)	$(LINKER) -o $@ $(OBJS) $(LIBRARIES) $(CUDA_LIBRARIES)

$(OBJDIR)quadtree.o				: quadtree.cpp
	$(VERBOSE)	$(CXX) $(CXXFLAGS) $(CUDA_INCLUDES) -o $@ -c quadtree.cpp


compile:	create_obj_dir $(OBJS)
create_obj_dir: 
	$(VERBOSE)	mkdir -p $(OBJDIR)
create_bin_dir:
	$(VERBOSE)	mkdir -p $(BINDIR)
	
clean:
	$(VERBOSE)	rm -vf $(OBJS)
tidy:	clean
	$(VERBOSE)	rm -vf $(TARGET)
clobber:	tidy
	$(VERBOSE)	rm -rvf obj/
	$(VERBOSE)	rm -rvf bin/
