
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
	cluster.o \
	new_aca.o \
	aca_b.o 

	
CC = gcc
CFLAGS = -Wall -Wextra -pedantic -g -DUSE_BLAS
LDLIBS = -lopenblas -lm -llapacke

 test: test.c $(OBJECTS)
	$(CC) $(CFLAGS)  test.c $(OBJECTS) -o $@ $(LDLIBS)

.PHONY: clean

clean:
	rm -f *.o  test
