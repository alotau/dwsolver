## Process this file with automake to produce Makefile.in ##

#INCLUDES = -I$(srcdir)/../include

# DW core files (everything except the CLI entry point and the new solver API).
DW_CORE_SOURCES = \
dw_blas.c \
dw_support.c \
dw_phases.c \
dw_subprob.c \
dw_rounding.c \
dw_globals.c

# -----------------------------------------------------------------------
# noinst convenience library: exposes dw_blas symbols for tests.
# Built without -fvisibility=hidden so test_blas can link against it
# without needing symbols exported from the full shared library.
# -----------------------------------------------------------------------
noinst_LTLIBRARIES = libdwblas_internal.la
libdwblas_internal_la_SOURCES = dw_blas.c
libdwblas_internal_la_CFLAGS  = $(AM_CFLAGS) -I$(srcdir)

# -----------------------------------------------------------------------
# Shared library: libdwsolver
# -----------------------------------------------------------------------
lib_LTLIBRARIES = libdwsolver.la

libdwsolver_la_SOURCES  = $(DW_CORE_SOURCES) dw_solver.c
libdwsolver_la_CPPFLAGS = $(AM_CPPFLAGS) -DDWSOLVER_BUILDING_LIB $(GLPK_CFLAGS)
libdwsolver_la_LIBADD   = $(GLPK_LIBS)
libdwsolver_la_CFLAGS   = $(AM_CFLAGS) -fvisibility=hidden
libdwsolver_la_LDFLAGS  = -version-info 0:0:0

# -----------------------------------------------------------------------
# CLI binary: dwsolver
# -----------------------------------------------------------------------
bin_PROGRAMS = dwsolver

dwsolver_SOURCES  = dw_main.c
dwsolver_CPPFLAGS = $(AM_CPPFLAGS) -DDWSOLVER_STATIC
dwsolver_CFLAGS   = $(AM_CFLAGS)
dwsolver_LDADD    = libdwsolver.la

# -----------------------------------------------------------------------
# Installed public header
# -----------------------------------------------------------------------
nobase_include_HEADERS = dw_solver.h

EXTRA_DIST = \
dw_blas.h \
dw.h \
dw_phases.h \
dw_rounding.h \
dw_subprob.h \
dw_support.h

## eof ##
