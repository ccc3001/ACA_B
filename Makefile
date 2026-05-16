
#
# Stand-alone Makefile
#

#export OPENBLAS_NUM_THREADS=1
#export OMP_NUM_THREADS=1


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
	aca_b.o \
	kernel_functions.o


CC = gcc
CFLAGS = -Wall -Wextra -pedantic -g -DUSE_BLAS
LDLIBS = -lopenblas -lm 
CFLAGS += $(shell pkg-config --cflags openblas)
LDFLAGS += $(shell pkg-config --libs openblas)
 test: test.c $(OBJECTS)
	$(CC) $(CFLAGS)  test.c $(OBJECTS) -o $@ $(LDLIBS)

.PHONY: clean

clean:
	rm -f *.o  test
