#!/bin/bash
#
# Set the following variables to the proper paths. To get help with the correct locations
# issue the following commands:
#   python-config --includes
#   python-config --ldflags
#
# PYTHON_LNAME is the name of the python library in the $PYTHON_LIB directory. If the
# library name is "libpython3.8", the variable should be set to "python3.8".
#
# It is assumed that you already built the LibGAL library. 
#
# If everything is successful you should see the "psgp4gm.so" library in this directory.
# Note, that it is OKay to use the ".so" extension on both Mac OS and Linux.  
#
PYTHON_INCLUDE="/opt/local/Library/Frameworks/Python.framework/Versions/3.8/include/python3.8"
PYTHON_LIB="/opt/local/Library/Frameworks/Python.framework/Versions/3.8/lib/python3.8/config-3.8-darwin"
PYTHON_LNAME="python3.8"
LIBGAL_INCLUDE="/opt/libgal-0.6.0/include"
LIBGAL_LIB="/opt/libgal-0.6.0/lib"
#
#
#
echo $PYTHON_INCLUDE
cython psgp4gm.pyx
gcc -c -fPIC -I$PYTHON_INCLUDE -I$LIBGAL_INCLUDE psgp4gm.c
gcc -shared psgp4gm.o -lgal -l$PYTHON_LNAME -L$LIBGAL_LIB -L$PYTHON_LIB -o psgp4gm.so
