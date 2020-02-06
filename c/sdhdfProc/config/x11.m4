# LIB_X11([search paths, [ACTION-IF-FOUND [,ACTION-IF-NOT-FOUND]]])
# ----------------------------------------------------------
AC_DEFUN([LIB_X11],
[
  AC_PROVIDE([LIB_X11])
  AC_REQUIRE([AC_PATH_XTRA])

  AC_MSG_CHECKING([for libX11.so])

  x11_found="no"

  if test "x$X_LIBS" != "x"; then
    x11_found="yes"
  else
    lib_x11_path_list="$1 /opt/X11/lib /usr/X11R6/lib"
    for dir in $lib_x11_path_list; do
      if test -r "$dir/libX11.so"; then
        X_LIBS="-L$dir"
        x11_found="yes"
      fi
    done
  fi

  AC_MSG_RESULT($x11_found)
  AC_MSG_RESULT([X_LIBS=$X_LIBS])
  AC_MSG_RESULT([X_PRE_LIBS=$X_PRE_LIBS])
  AC_MSG_RESULT([X_EXTRA_LIBS=$X_EXTRA_LIBS])
  AC_MSG_RESULT([X_CFLAGS=$X_CFLAGS])


])
