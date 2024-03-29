dnl Process this file with autoconf to produce a configure script

AC_INIT([dwsolver], [1.2.1], [Joseph.L.Rios@nasa.gov])

dnl AC_CONFIG_SRCDIR([include/glpk.h])

AC_CONFIG_MACRO_DIR([m4])

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

AC_ARG_WITH(gmp,
AC_HELP_STRING([--with-gmp],
   [use GNU MP bignum library [[default=no]]]),
   [case $withval in
      yes | no) ;;
      *) AC_MSG_ERROR([invalid value `$withval' for --with-gmp]);;
      esac],
   [with_gmp=no])

AC_ARG_WITH(zlib,
AC_HELP_STRING([--with-zlib],
   [use zlib data compression library [[default=no]]]),
   [case $withval in
      yes | no) ;;
      *) AC_MSG_ERROR([invalid value `$withval' for --with-zlib]);;
      esac],
   [with_zlib=no])

AC_ARG_ENABLE(dl,
AC_HELP_STRING([--enable-dl],
   [enable shared library support [[default=no]]]),
   [case $enableval in
      yes | ltdl | dlfcn | no) ;;
      *) AC_MSG_ERROR([invalid value `$enableval' for --enable-dl]);;
      esac],
   [enable_dl=no])

AC_ARG_ENABLE(odbc,
AC_HELP_STRING([--enable-odbc],
   [enable MathProg ODBC support [[default=no]]]),
   [case $enableval in
      yes | unix | no) ;;
      *) AC_MSG_ERROR([invalid value `$enableval' for --enable-odbc]);;
      esac],
   [enable_odbc=no])

AC_ARG_ENABLE(mysql,
AC_HELP_STRING([--enable-mysql],
   [enable MathProg MySQL support [[default=no]]]),
   [case $enableval in
      yes | no) ;;
      *) AC_MSG_ERROR([invalid value `$enableval' for --enable-mysql]);;
      esac],
   [enable_mysql=no])

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

AC_MSG_CHECKING([whether to use GNU MP bignum library])
if test "$with_gmp" = "yes"; then
   AC_MSG_RESULT([yes])
   AC_DEFINE([HAVE_GMP], [1], [N/A])
   LIBS="-lgmp $LIBS"
else
   AC_MSG_RESULT([no])
fi

AC_MSG_CHECKING([whether to use zlib data compression library])
if test "$with_zlib" = "yes"; then
   AC_MSG_RESULT([yes])
   AC_DEFINE([HAVE_ZLIB], [1], [N/A])
   LIBS="-lz $LIBS"
else
   AC_MSG_RESULT([no])
fi

AC_MSG_CHECKING([whether to enable shared library support])
if test "$enable_dl" = "yes"; then
   AC_MSG_RESULT([ltdl])
   AC_DEFINE([HAVE_LTDL], [1], [N/A])
   LIBS="-lltdl $LIBS"
elif test "$enable_dl" = "ltdl"; then
   AC_MSG_RESULT([ltdl])
   AC_DEFINE([HAVE_LTDL], [1], [N/A])
   LIBS="-lltdl $LIBS"
elif test "$enable_dl" = "dlfcn"; then
   AC_MSG_RESULT([dlfcn])
   AC_DEFINE([HAVE_DLFCN], [1], [N/A])
else
   AC_MSG_RESULT([no])
fi

case $host_os in
   darwin* | macosx*)
      LIBIODBC="libiodbc.dylib"
      LIBODBC="libodbc.dylib"
      LIBMYSQL="libmysqlclient.dylib"
      CFLAGS="-DUSE_NAMED_SEMAPHORES $CFLAGS"
      ;;
   *)
      LIBIODBC="libiodbc.so"
      LIBODBC="libodbc.so"
      LIBMYSQL="libmysqlclient.so"
      ;;
esac

AC_MSG_CHECKING([whether to enable MathProg ODBC support])
if test "$enable_odbc" = "yes"; then
   if test "$enable_dl" = "no"; then
      AC_MSG_ERROR([--enable-odbc requires --enable-dl])
   fi
   AC_MSG_RESULT([yes])
   AC_DEFINE_UNQUOTED([ODBC_DLNAME], ["$LIBIODBC"], [N/A])
elif test "$enable_odbc" = "unix"; then
   if test "$enable_dl" = "no"; then
      AC_MSG_ERROR([--enable-odbc requires --enable-dl])
   fi
   AC_MSG_RESULT([unix])
   AC_DEFINE_UNQUOTED([ODBC_DLNAME], ["$LIBODBC"], [N/A])
else
   AC_MSG_RESULT([no])
fi

AC_MSG_CHECKING([whether to enable MathProg MySQL support])
if test "$enable_mysql" = "yes"; then
   if test "$enable_dl" = "no"; then
      AC_MSG_ERROR([--enable-mysql requires --enable-dl])
   fi
   AC_MSG_RESULT([yes])
   CPPFLAGS="-I/usr/include/mysql $CPPFLAGS"
   AC_DEFINE_UNQUOTED([MYSQL_DLNAME], ["$LIBMYSQL"], [N/A])
else
   AC_MSG_RESULT([no])
fi

AC_CONFIG_FILES(
dnl   [include/Makefile src/Makefile Makefile])
	[src/Makefile Makefile])
AC_OUTPUT

dnl eof
