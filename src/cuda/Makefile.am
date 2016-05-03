#libPFScan_ProcessGPU.cu.o:
#       nvcc -std=c++11 -m64 -g PFScan_ProcessGPU.cu  PFScan.h $(CPPFLAGS) $(LDFLAGS)
#libcudafiles_la_LINK= $(LIBTOOL) --mode=link $(CXX) -o $@ $(CUDA_LDFLAGS) $(CUDA_LIBS)

noinst_LIBRARIES = PFScan_ProcessGPU.a
#noinst_LIBRARIES = libPFScan_ProcessGPU.a
#noinst_LTLIBRARIES = PFScan_ProcessGPU.la

#CC = nvcc
#AM_CXXFLAGS = -std=c++11 -m64 -g
#AM_CPPFLAGS = -I$(abs_top_builddir)/include

#noinst_LIBRARIES = libPFScan_ProcessGPU.a
PFScan_ProcessGPU_a_SOURCES = PFScan.h PFScan_ProcessGPU.cu  
PFScan_ProcessGPU_a_CFLAGS = -I$(abs_top_builddir)/include -I${ROOTSYS}/include/root -I${CUDAHOME}/include
PFScan_ProcessGPU_a_AR = nvcc -std=c++11 -arch=sm_20 -m64 -g -lib -o

.cu.o: PFScan.h
	nvcc $(PFScan_ProcessGPU_a_CFLAGS) -std=c++11 -arch=sm_20 -m64 -g -dc -o $@ $<

#libPFScan_ProcessGPU_la_LIBTOOLFLAGS = --mode=compile nvcc -std=c++11 -m64 -g -I$(abs_top_builddir)/include
#libPFScan_ProcessGPU_la_CC = nvcc
#libPFScan_ProcessGPU_la_CFLAGS = -std=c++11 -m64 -g
#libPFScan_ProcessGPU_la_CPPFLAGS = -I$(abs_top_builddir)/include

