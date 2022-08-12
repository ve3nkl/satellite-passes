#!/bin/bash
#
# Check that the following variables are set correctly for your environment.
#
# PYTHON_LNAME is the name of the python library in the $PYTHON_LIB directory. If the
# library name is "libpython3.9", the variable should be set to "python3.9".
#
# It is assumed that you already built the LibGAL library. 
#
# If everything is successful you should see the "pygal.so" library in this directory.
# Note, that it is OKay to use the ".so" extension on both Mac OS and Linux.  
#
PYTHON_LNAME="python3.9"
LIBGAL_INCLUDE="/opt/libgal-0.6.0/include"
LIBGAL_LIB="/opt/libgal-0.6.0/lib"
#
PYTHON_INCLUDE=$(python-config --includes)
PYTHON_LIB=$(python-config --ldflags)
#
#
echo "Calling cython ..."
cython pygal.pyx
if [ $? -eq 0 ]
then
  echo "Compiling pygal.c ..."
  gcc -c -fPIC $PYTHON_INCLUDE -I$LIBGAL_INCLUDE pygal.c
  if [ $? -eq 0 ]
  then
    echo "Building pygal.so ..."
    gcc -shared pygal.o -lgal -l$PYTHON_LNAME -L$LIBGAL_LIB $PYTHON_LIB -o pygal.so
    if [ $? -eq 0 ]
    then 
      echo "... built successfully."
    fi
  fi
fi
