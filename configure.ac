dnl Process this file with autoconf to produce a configure script.

dnl We need at least autoconf 2.69 for this configure.ac to work.
AC_PREREQ([2.69])
AC_INIT([cddplus], [0.77], [julian.rueth@fsfe.org])
AM_INIT_AUTOMAKE([-Wall -Werror foreign])

AC_PROG_CXX
AC_CANONICAL_HOST
AC_PROG_INSTALL
LT_INIT

# We compile against GMP.
AC_CHECK_HEADER(gmp.h, , [AC_MSG_ERROR([GMP header not found])])
AC_CHECK_LIB(gmp, __gmpz_add, [], [AC_MSG_ERROR([libgmp not found])])

AC_CONFIG_FILES([Makefile])

AC_OUTPUT
