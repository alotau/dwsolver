dnl Process this file with autoconf to produce a configure script

AC_INIT([dwsolver], [1.2.1], [Joseph.L.Rios@nasa.gov])

dnl AC_CONFIG_SRCDIR([include/glpk.h])

AC_CONFIG_MACRO_DIR([m4])

AC_CONFIG_AUX_DIR([build-aux])
AM_INIT_AUTOMAKE

AC_CONFIG_HEADERS([config.h])

AC_ARG_ENABLE(named_semaphores,
AC_HELP_STRING([--enable-named-semaphores],
   [enable named semaphores (probably needed on Mac) [[default=no]]]),
   [case $enableval in
      yes | no) ;;
      *) AC_MSG_ERROR([invalid value `$withval' for --enable-named-semaphores]);;
      esac],
   [enable_named_semaphores=no])
   
AC_ARG_ENABLE(recursive_mutex,
AC_HELP_STRING([--enable-recursive-mutex],
   [enable use of recursive mutices (probably needed on Mac) [[default=no]]]),
   [case $enableval in
      yes | no) ;;
      *) AC_MSG_ERROR([invalid value `$withval' for --enable-recursive-mutex]);;
      esac],
   [enable_recursive_mutex=no])

dnl Disable unnecessary libtool tests
define([AC_LIBTOOL_LANG_CXX_CONFIG], [:])
define([AC_LIBTOOL_LANG_F77_CONFIG], [:])
define([AC_LIBTOOL_LANG_GCJ_CONFIG], [:])

dnl Check for programs
AC_PROG_CC
AC_PROG_INSTALL
AC_PROG_LIBTOOL
AM_PROG_CC_C_O

dnl Check for math library
AC_CHECK_LIB([m], [exp])
AC_CHECK_LIB([pthread], [sem_wait])

dnl Require system GLPK >= 4.65
PKG_CHECK_MODULES([GLPK], [glpk >= 4.65],
  [],
  [AC_MSG_ERROR([GLPK >= 4.65 is required. Install libglpk-dev (Debian/Ubuntu), glpk-devel (Fedora/RHEL), or brew install glpk (macOS).])])

dnl Check for <sys/time.h> header
AC_CHECK_HEADER([sys/time.h],
   AC_DEFINE([HAVE_SYS_TIME_H], [1], [N/A]))

dnl Check for gettimeofday function
AC_CHECK_FUNC([gettimeofday],
   AC_DEFINE([HAVE_GETTIMEOFDAY], [1], [N/A]))
   
dnl Check if we should use named semaphores.  Can do this based on OS name?
AC_MSG_CHECKING([whether to use named semaphores])
if test "$enable_named_semaphores" = "yes"; then
   AC_MSG_RESULT([yes])
   CFLAGS="-DUSE_NAMED_SEMAPHORES $CFLAGS"
else
   AC_MSG_RESULT([no])
fi

AC_MSG_CHECKING([whether to use recursive mutices])
if test "$enable_recursive_mutex" = "yes"; then
   AC_MSG_RESULT([yes])
   CFLAGS="-DUSE_RECURSIVE_MUTEX $CFLAGS"
else
   AC_MSG_RESULT([no])
fi

AC_CONFIG_FILES(
dnl   [include/Makefile src/Makefile Makefile])
	[src/Makefile Makefile tests/Makefile dwsolver.pc])
AC_OUTPUT

dnl eof
