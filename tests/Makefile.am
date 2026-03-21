## Process this file with automake to produce Makefile.in ##

TESTS           = test_blas
check_PROGRAMS  = test_blas

test_blas_SOURCES = test_blas.c
test_blas_CFLAGS  = -I$(top_srcdir)/src
test_blas_LDADD   = $(top_builddir)/src/dw_blas.o

## eof ##
