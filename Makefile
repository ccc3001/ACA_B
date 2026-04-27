
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

example_11: example_11.c $(OBJECTS)
	$(CC) $(CFLAGS) example_11.c $(OBJECTS) -o $@ $(LDLIBS)

.PHONY: clean

clean:
	rm -f *.o example_11
