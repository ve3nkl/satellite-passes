#!/bin/bash

PYTHON_INCLUDES=$(python-config --includes)
PYTHON_LDFLAGS=$(python-config --ldflags)
cython psgp4gm.pyx
gcc -c -fPIC $PYTHON_INCLUDES -I/opt/libgal-0.6.0/include psgp4gm.c
gcc -shared psgp4gm.o -lgal -L/opt/libgal-0.6.0/lib -lpython3.8 $PYTHON_LDFLAGS -o psgp4gm.so
