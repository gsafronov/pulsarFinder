TOP = $(shell pwd)
TARGET = run_dmscan

CC = g++
CU = nvcc

CCFLAGS = -std=c++11 -pthread -m64 -g
CUFLAGS = -std=c++11 -m64 -g

SRCDIR = src
CPPOBJECTS = $(SRCDIR)/DMScan_cppfunc.cc src/DMAnalysis_cppfunc.cc
CUOBJECTS = $(SRCDIR)/DMScan_cudafunc.cu

LDIR = ./lib
INCLUDEDIR = -I$(ROOTSYS)/include/root -I. -I./include -I$(CUDAHOME)/include
LIBDIR = -L$(ROOTSYS)/lib/root -L. -L./lib

LDFLAGS = -lCore \
        -lRIO \
        -lNet \
        -lHist \
        -lGraf \
        -lGraf3d \
        -lGpad \
        -lTree \
        -lRint \
        -lPostscript \
        -lMatrix \
        -lPhysics \
        -lMathCore \
        -lThread \
        -lz \
        -lGui \
        -pthread \
        -lm \
        -ldl \
	-rdynamic \
	-lstdc++

LIBCUDA = -L$(CUDAHOME)/lib64 -L. -L./lib

LDFLAGS_CUDA = -lcudart -lcufft
	

all: run_dmscan run_dman run_readRAW

run_readRAW: readRAW
	$(CC) $(CCFLAGS) $(LDIR)/readRAW.o -o bin/readRAW $(INCLUDEDIR) $(LIBDIR) $(LIBCUDA) $(LDFLAGS)  

readRAW: 
	$(CC) $(CCFLAGS) src/readRAW.cc -c -o $(LDIR)/readRAW.o $(INCLUDEDIR) $(LIBDIR) $(LDFLAGS)

run_dman: dman DMScan_cppfunc DMAnalysis_cppfunc DMScan_cudafunc
	$(CC) $(CCFLAGS) $(LDIR)/dman.o $(LDIR)/DMScan_cppfunc.o $(LDIR)/DMAnalysis_cppfunc.o $(LDIR)/DMScan_cudafunc.o -o bin/run_dman $(INCLUDEDIR) $(LIBDIR) $(LIBCUDA) $(LDFLAGS) $(LDFLAGS_CUDA)

run_dmscan: dmscan DMScan_cppfunc DMAnalysis_cppfunc DMScan_cudafunc
	$(CC) $(CCFLAGS) $(LDIR)/dmscan.o $(LDIR)/DMScan_cppfunc.o $(LDIR)/DMAnalysis_cppfunc.o $(LDIR)/DMScan_cudafunc.o -o bin/run_dmscan $(INCLUDEDIR) $(LIBDIR) $(LIBCUDA) $(LDFLAGS) $(LDFLAGS_CUDA)

dman: 
	$(CC) $(CCFLAGS) src/run_dman.cc -c -o $(LDIR)/dman.o $(INCLUDEDIR) $(LIBDIR) $(LIBCUDA) $(LDFLAGS) $(LDFLAGS_CUDA)

dmscan:  
	$(CC) $(CCFLAGS) src/run_dmscan.cc -c -o $(LDIR)/dmscan.o $(INCLUDEDIR) $(LIBDIR) $(LIBCUDA) $(LDFLAGS)

DMScan_cppfunc:
	$(CC) $(CCFLAGS) $(SRCDIR)/DMScan_cppfunc.cc -c -o $(LDIR)/DMScan_cppfunc.o $(INCLUDEDIR) $(LIBDIR) $(LDFLAGS)

DMAnalysis_cppfunc:
	$(CC) $(CCFLAGS) $(SRCDIR)/DMAnalysis_cppfunc.cc -c -o $(LDIR)/DMAnalysis_cppfunc.o $(INCLUDEDIR) $(LIBDIR) $(LDFLAGS)

DMScan_cudafunc:
	$(CU) $(CUFLAGS) $(SRCDIR)/DMScan_cudafunc.cu -c -o $(LDIR)/DMScan_cudafunc.o $(INCLUDEDIR) 

clean:
	rm -rf $(LDIR)/*o run_dmscan.exe
