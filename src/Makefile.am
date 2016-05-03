bin_PROGRAMS = readRAW runPFinder  
AM_CPPFLAGS = -I$(abs_top_builddir)/include
AM_LDFLAGS = -L$(abs_top_builddir)/src/cuda -L.
#-lPFScan_ProcessGPU
#SUBDIRS = cuda 
CC = g++
AM_CXXFLAGS = -std=c++11 -pthread -m64 -g
readRAW_SOURCES = readRAW.cc
SUBDIRS = cuda

runPFinder_SOURCES = runPFinder.cc PFScan_ProcessCPU.cc PFScan_InitClose.cc PFRun.cc PFAnalysis.cc PFScan.h PFRun.h PFAnalysis.h 
runPFinder_LDADD = pfscan_cudalink.o cuda/PFScan_ProcessGPU.a 
#pfscan_cudalink.o cuda/PFScan_ProcessGPU.a -lcuda

pfscan_cudalink.o: runPFinder.o cuda/PFScan_ProcessGPU.a
	nvcc -std=c++11 -arch=sm_20 -m64 -g -dlink runPFinder.o cuda/PFScan_ProcessGPU.a -o pfscan_cudalink.o
