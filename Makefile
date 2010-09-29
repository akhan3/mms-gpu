CXX         := 	g++
CC         	:=	gcc
LINKER     	:= 	g++ -fPIC

# CUDA compiler, includes and libraries
NVCC 				:= 	nvcc
CUDA_INSTALL_PATH 	:= 	/usr/local/cuda
CUDA_SDK_PATH 		:= 	$(HOME)/NVIDIA_GPU_Computing_SDK
CUDA_INCLUDES   	:= 	-I$(CUDA_INSTALL_PATH)/include -I$(CUDA_SDK_PATH)/C/common/inc
CUDA_LIBRARIES 		:= 	-L$(CUDA_INSTALL_PATH)/lib64 -lcuda -lcudart \
						-L$(CUDA_SDK_PATH)/C/lib/ -lcutil_x86_64

#INCPATH		:= -I./ $(MATLAB_INCLUDES)
#LIBPATH		:= $(MATLAB_LIBRARIES)
DEFINES		:=
COMMONFLAGS	:=
LIBRARIES	:=

ifeq ($(use_freeimage),1)
	LIBRARIES	+= -lfreeimage
	DEFINES 	+= -DUSE_FREEIMAGE
endif

COMMONFLAGS		+= $(DEFINES)
ifeq ($(dbg),1)
	COMMONFLAGS += -g
	NVCCFLAGS	+= -g -G
else
	COMMONFLAGS += -O3
	NVCCFLAGS	+= -O3
endif

ifeq ($(prof),1)
	COMMONFLAGS += -pg
endif
ifeq ($(omp),0)
	COMMONFLAGS +=
else
	COMMONFLAGS += -fopenmp
endif

NVCCFLAGS   += 	--compiler-options "-fno-strict-aliasing $(COMMONFLAGS) -W" $(CUDA_INCLUDES) \
						-arch=sm_20
						#-gencode=arch=compute_20,code=sm_20 \
						#-gencode=arch=compute_20,code=compute_20

CXXFLAGS    :=	-Wall -W $(INCPATH) $(COMMONFLAGS)
OBJS        :=	fmm_kernel.o \
				helper_functions.o \
				vector_functions.o \
				ode_functions.o \
				main.o

TARGET    	:= main

all: $(TARGET)

$(TARGET):	$(OBJS)
	@echo "******** Linking ********"
	$(LINKER) -o $@ $(OBJS) $(COMMONFLAGS) $(LIBRARIES) $(CUDA_LIBRARIES)

fmm_kernel.o	: fmm_kernel.cu fmm_calc.cu potential_calc.cu Box.* Queue.* Cmpx.* Vector3.* numerics.*  mydefs.hpp
	$(NVCC) $(NVCCFLAGS) -o $@ -c $<
helper_functions.o 	: helper_functions.cpp helper_functions.hpp mydefs.hpp
vector_functions.o 	: vector_functions.cpp vector_functions.hpp mydefs.hpp
ode_functions.o 	: ode_functions.cpp ode_functions.hpp mydefs.hpp
main.o				: main.cpp mydefs.hpp

clean:
	rm -f $(OBJS) gmon.out $(TARGET)
