#libDMScan_cudafunc.cu.o:
#       nvcc -std=c++11 -m64 -g DMScan_cudafunc.cu  DMScan.h $(CPPFLAGS) $(LDFLAGS)
#libcudafiles_la_LINK= $(LIBTOOL) --mode=link $(CXX) -o $@ $(CUDA_LDFLAGS) $(CUDA_LIBS)

noinst_LIBRARIES = DMScan_cudafunc.a
#noinst_LIBRARIES = libDMScan_cudafunc.a
#noinst_LTLIBRARIES = DMScan_cudafunc.la

#CC = nvcc
#AM_CXXFLAGS = -std=c++11 -m64 -g
#AM_CPPFLAGS = -I$(abs_top_builddir)/include

#noinst_LIBRARIES = libDMScan_cudafunc.a
DMScan_cudafunc_a_SOURCES = DMScan.h DMScan_cudafunc.cu  
DMScan_cudafunc_a_CFLAGS = -I$(abs_top_builddir)/include -I${ROOTSYS}/include/root -I${CUDAHOME}/include
DMScan_cudafunc_a_AR = nvcc -std=c++11 -arch=sm_20 -m64 -g -lib -o

.cu.o: DMScan.h
	nvcc $(DMScan_cudafunc_a_CFLAGS) -std=c++11 -arch=sm_20 -m64 -g -dc -o $@ $<

#libDMScan_cudafunc_la_LIBTOOLFLAGS = --mode=compile nvcc -std=c++11 -m64 -g -I$(abs_top_builddir)/include
#libDMScan_cudafunc_la_CC = nvcc
#libDMScan_cudafunc_la_CFLAGS = -std=c++11 -m64 -g
#libDMScan_cudafunc_la_CPPFLAGS = -I$(abs_top_builddir)/include

