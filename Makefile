
#
# Stand-alone Makefile
#

OBJECTS = \
	aca.o \
	basic.o \
	blas.o \
	fullmatrix.o \
	rkmatrix.o \
	supermatrix.o \
	interpolation.o \
	cluster.o
	
CC = gcc
CFLAGS = -Wall -Wextra -pedantic -g -DUSE_BLAS
LDLIBS = -lopenblas -lm

 kernel_functions: kernel_functions.c $(OBJECTS)
	$(CC) $(CFLAGS)  kernel_functions.c $(OBJECTS) -o $@ $(LDLIBS)

.PHONY: clean

clean:
	rm -f *.o  kernel_functions
