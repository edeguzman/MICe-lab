# Require autoconf 2.57 or newer.
AC_PREREQ([2.57])
AC_INIT([MICe-lab],[0.14],[Matthijs van Eede matthijs.vaneede@sickkids.ca])
#foreign tells Automake that this package does not follow GNU standards.
AM_INIT_AUTOMAKE([-Wall -Werror foreign])
AC_PROG_CC
AC_PROG_CXX

PKG_CHECK_EXISTS(opencv, [
  CV_LIB_EXTRA="$(pkg-config --libs-only-L opencv | cut -f1 -d' ' | sed 's|-L||' | sed 's|lib/?$|share/OpenCV/3rdparty/lib|')"
  LDFLAGS="-L$CV_LIB_EXTRA"
  CPPFLAGS="$(pkg-config --cflags --libs opencv)"],
  [AC_MSG_ERROR([OpenCV not found.])]
)

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

