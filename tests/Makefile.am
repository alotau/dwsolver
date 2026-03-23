## Process this file with automake to produce Makefile.in ##

TESTS           = test_blas test_lib_api
check_PROGRAMS  = test_blas test_lib_api

test_blas_SOURCES = test_blas.c
test_blas_CFLAGS  = -I$(top_srcdir)/src
test_blas_LDADD   = $(top_builddir)/src/dw_blas.o

test_lib_api_SOURCES  = test_lib_api.c
test_lib_api_CPPFLAGS = -I$(top_srcdir)/src -DEXAMPLE_DIR='"$(top_srcdir)/examples"' -DDWSOLVER_STATIC
test_lib_api_LDADD    = $(top_builddir)/src/libdwsolver.la

EXTRA_DIST = \
	dw-tests.sh \
	test_guidefile.sh \
	fixtures

# Files written by dw-tests.sh into the build tree during 'make check'.
# Listed here so 'make distclean' (and hence 'make distcheck') leaves the
# tree clean enough to pass the distcleancheck.
DISTCLEANFILES = \
	relaxed_solution \
	rs_sorted \
	out_obj.txt \
	src/dw_blas.Po

## eof ##
