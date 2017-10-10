# Makefile
# Generic Makefile for making cuda programs
#
BIN					:= mdis_trunk
# flags
CUDA_INSTALL_PATH	:= /usr/local/cuda
OBJDIR				:= obj
INCLUDES			+= -I$(CUDA_INSTALL_PATH)/include -I.
LIBS				:= -L$(CUDA_INSTALL_PATH)/lib64
CFLAGS				:= -O0 -g
LDFLAGS				:= -lrt -lm -lcudart
# compilers
#NVCC				:= $(CUDA_INSTALL_PATH)/bin/nvcc --compiler-options -fpermissive -arch sm_20 --ptxas-options=-v
NVCC				:= $(CUDA_INSTALL_PATH)/bin/nvcc --compiler-options -fpermissive --ptxas-options=-v -use_fast_math 
CC					:= g++
LINKER				:= g++ -fPIC -g
# files
CPP_SOURCES			:= \
  Core/main.cpp \
  Core/parameters.cpp \
  IO/configreader.cpp \
  IO/pdbio.cpp \
  IO/xyzio.cpp \
  IO/psfio.cpp \
  IO/topio.cpp \
  IO/dcdio.cpp \
  IO/quadratureio.cpp \
  IO/gbswpbradiiio.cpp \
  IO/ffreader.cpp \
  Core/ffutils.cpp \
  Core/topology.cpp \
  Core/atomTypes.cpp \
  Core/angleTypes.cpp \
  Util/mdunitsconverter.cpp \
  Potentials/sasa.cpp \
  Util/timer.cpp \
  IO/filetypes.cpp \
  Util/memory.cpp \
  Util/mystl.cpp \
  Util/wrapper.cpp \
  Util/atomfilter.cpp \
  ConstrAlgorithms/quern/quern.cpp \
  ConstrAlgorithms/quern/quern_order.cpp \
  ConstrAlgorithms/quern/quern_solver.cpp \
  ConstrAlgorithms/quern/quern_factorization.cpp

CU_SOURCES			:= Core/md.cu
CPP_OBJS				:= $(patsubst %.cpp, $(OBJDIR)/%.cpp.o, $(CPP_SOURCES))
CU_OBJS				:= $(patsubst %.cu, $(OBJDIR)/%.cu.o, $(CU_SOURCES))
 
$(BIN): makedirs clean $(CU_OBJS)
	$(LINKER) -o $(BIN) $(CU_OBJS) $(CPP_SOURCES) $(LDFLAGS) $(INCLUDES) $(LIBS)
 
$(OBJDIR)/%.c.o: $(CPP_SOURCES)
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $<
 
$(OBJDIR)/%.cu.o: $(CU_SOURCES)
	$(NVCC) $(INCLUDES) -o $@ -c $<
 
makedirs:
	mkdir -p $(OBJDIR)/Util
	mkdir -p $(OBJDIR)/Core
	mkdir -p $(OBJDIR)/IO
	mkdir -p $(OBJDIR)/Analysis

run: $(BIN)
	LD_LIBRARY_PATH=$(CUDA_INSTALL_PATH)/lib ./$(BIN)
 
clean:
	rm -f $(BIN) $(OBJDIR)/Util/*.o
	rm -f $(BIN) $(OBJDIR)/Core/*.o
	rm -f $(BIN) $(OBJDIR)/IO/*.o
	
install:
	cp $(BIN) /usr/bin/$(BIN)
	
