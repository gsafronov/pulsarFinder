#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([pulsarFinder], [1.1.1])
AM_INIT_AUTOMAKE([foreign])
AC_CONFIG_SRCDIR([src/PFScan_ProcessCPU.cc])
#AC_CONFIG_HEADERS([header.h])

#check if ROOTSYS is set:
if [ test -n "$ROOTSYS" && test -n "$CUDAHOME" ]
then
    LDFLAGS="$LDFLAGS -L${ROOTSYS}/lib/root -L${CUDAHOME}/lib64 -L${CUDAHOME}/lib -L."
    CPPFLAGS="-I${abs_top_builddir}/include -I${ROOTSYS}/include/root -I${CUDAHOME}/include"
elif [ ! test -n "$ROOTSYS" ]
then
    echo 'ROOTSYS variable is not set. If in bash type export ROOTSYS=<head directory of root>, in csh type setenv ROOTSYS <head directory of root>. Exiting..'
    exit
elif [ ! test -n "$CUDAHOME" ]
then
    echo 'CUDAHOME variable is not set. If in bash type export ROOTSYS=<head directory of root>, in csh type setenv ROOTSYS <head directory of root>. Exiting..'
    exit
fi

#AC_ENABLE_SHARED
#AC_DISABLE_STATIC

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC
AC_PROG_RANLIB
AM_PROG_AR
#LT_INIT
#AC_PROG_LIBTOOL
AC_CHECK_PROG( NVCCYES, nvcc, true, false, [${CUDAHOME}/bin])

if [ ! $NVCCYES ] 
then
    echo 'nvcc was not found in $CUDAHOME/bin'
    exit
fi

LDFLAGS="$LDFLAGS -L${ROOTSYS}/lib/root -L${CUDAHOME}/lib64"

# Checks for libraries.
# FIXME: Replace `main' with a function in `-lCore':
AC_CHECK_LIB([Core], [main])
# FIXME: Replace `main' with a function in `-lGpad':
AC_CHECK_LIB([Gpad], [main])
# FIXME: Replace `main' with a function in `-lGraf':
AC_CHECK_LIB([Graf], [main])
# FIXME: Replace `main' with a function in `-lGraf3d':
AC_CHECK_LIB([Graf3d], [main])
# FIXME: Replace `main' with a function in `-lGui':
AC_CHECK_LIB([Gui], [main])
# FIXME: Replace `main' with a function in `-lHist':
AC_CHECK_LIB([Hist], [main])
# FIXME: Replace `main' with a function in `-lMathCore':
AC_CHECK_LIB([MathCore], [main])
# FIXME: Replace `main' with a function in `-lMatrix':
AC_CHECK_LIB([Matrix], [main])
# FIXME: Replace `main' with a function in `-lNet':
AC_CHECK_LIB([Net], [main])
# FIXME: Replace `main' with a function in `-lPhysics':
AC_CHECK_LIB([Physics], [main])
# FIXME: Replace `main' with a function in `-lPostscript':
AC_CHECK_LIB([Postscript], [main])
# FIXME: Replace `main' with a function in `-lRIO':
AC_CHECK_LIB([RIO], [main])
# FIXME: Replace `main' with a function in `-lRint':
AC_CHECK_LIB([Rint], [main])
# FIXME: Replace `main' with a function in `-lThread':
AC_CHECK_LIB([Thread], [main])
# FIXME: Replace `main' with a function in `-lTree':
AC_CHECK_LIB([Tree], [main])
# FIXME: Replace `main' with a function in `-lcudart':
AC_CHECK_LIB([cudart], [main])
# FIXME: Replace `main' with a function in `-lcufft':
AC_CHECK_LIB([cufft], [main])
# FIXME: Replace `main' with a function in `-ldl':
AC_CHECK_LIB([dl], [main])
# FIXME: Replace `main' with a function in `-lm':
AC_CHECK_LIB([m], [main])
# FIXME: Replace `main' with a function in `-lz':
AC_CHECK_LIB([z], [main])

# Checks for header files.
AC_CHECK_HEADERS([stdint.h stdlib.h string.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_CHECK_HEADER_STDBOOL
AC_TYPE_UINT32_T

# Checks for library functions.
AC_CHECK_FUNCS([floor mkdir pow sqrt])

AC_CONFIG_FILES([Makefile src/Makefile src/cuda/Makefile])
AC_OUTPUT

