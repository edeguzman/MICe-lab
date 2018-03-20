# Require autoconf 2.57 or newer.
AC_PREREQ([2.57])

AC_INIT([MICe-lab],[0.14],[Matthijs van Eede matthijs.vaneede@sickkids.ca])

#foreign tells Automake that this package does not follow GNU standards.
AM_INIT_AUTOMAKE([-Wall -Werror foreign])

AC_PROG_CC
AC_PROG_CXX

# The argument is the name of the generated header file.
# It is recommended to leave it as "config.h".
AC_CONFIG_HEADERS([config.h])

AC_CONFIG_FILES([
 Makefile
 src/Makefile
 perl/Makefile
])
# Checks for programs.
AC_PROG_INSTALL

# Finish up by writing output files.
AC_OUTPUT

