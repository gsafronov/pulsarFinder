bin_PROGRAMS = readRAW run_dmscan run_dman 
AM_CPPFLAGS = -I$(abs_top_builddir)/include
AM_LDFLAGS = -L$(abs_top_builddir)/src/cuda -L.
#-lDMScan_cudafunc
#SUBDIRS = cuda 
CC = g++
AM_CXXFLAGS = -std=c++11 -pthread -m64 -g
readRAW_SOURCES = readRAW.cc
SUBDIRS = cuda

run_dmscan_SOURCES = run_dmscan.cc DMScan_cppfunc.cc DMScan.h 
run_dmscan_LDADD = dmscan_cudalink.o cuda/DMScan_cudafunc.a 
#dmscan_cudalink.o cuda/DMScan_cudafunc.a -lcuda

run_dman_SOURCES = run_dman.cc DMScan_cppfunc.cc DMAnalysis_cppfunc.cc DMAnalysis.h DMScan.h
run_dman_LDADD = dman_cudalink.o cuda/DMScan_cudafunc.a 
#dman_cudalink.o cuda/DMScan_cudafunc.a -lcuda

dmscan_cudalink.o: run_dmscan.o cuda/DMScan_cudafunc.a
	nvcc -std=c++11 -arch=sm_20 -m64 -g -dlink run_dmscan.o cuda/DMScan_cudafunc.a -o dmscan_cudalink.o

dman_cudalink.o: run_dman.o cuda/DMScan_cudafunc.a
	nvcc -std=c++11 -arch=sm_20 -m64 -g -dlink run_dman.o cuda/DMScan_cudafunc.a -o dman_cudalink.o