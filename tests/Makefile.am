## Process this file with automake to produce Makefile.in ##

TESTS           = test_blas test_lib_api
check_PROGRAMS  = test_blas test_lib_api

test_blas_SOURCES = test_blas.c
test_blas_CFLAGS  = -I$(top_srcdir)/src
test_blas_LDADD   = $(top_builddir)/src/dw_blas.o

test_lib_api_SOURCES  = test_lib_api.c
test_lib_api_CPPFLAGS = -I$(top_srcdir)/src -DEXAMPLE_DIR='"$(top_srcdir)/examples"'
test_lib_api_LDADD    = $(top_builddir)/src/libdwsolver.la

## eof ##
