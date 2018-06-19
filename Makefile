

COMPILER = gcc
LIBRARIES = /Users/danafaiez/clapack/CLAPACK-3.2.1/tmglib_LINUX.a /Users/danafaiez/clapack/CLAPACK-3.2.1/lapack_LINUX.a /Users/danafaiez/clapack/CLAPACK-3.2.1/blas_LINUX.a /Users/danafaiez/clapack/CLAPACK-3.2.1/F2CLIBS/libf2c.a /usr/local/lib/libgsl.a /usr/local/lib/libgslcblas.a -lm -lc 
CFLAGS = -g -fno-standalone-debug #-O6# 
C_OBJECTS = $(patsubst %.c,%.o,$(wildcard *.c))
LIBDIR = 
#EFENCE =   -lefence

all: diag

ft: ft

$(C_OBJECTS): %.o: %.c
	$(COMPILER) -c $(CFLAGS)  $< -o $@


diag:  $(C_OBJECTS)
	$(COMPILER) $(CFLAGS) -o $@  $^  $(LIBDIR) $(LIBRARIES) $(EFENCE)

ft:  ft.o
	$(COMPILER) $(CFLAGS) -o $@  $^ -lm

det: determinant.o zgetrf_dbl_cmplx.o
	$(COMPILER) $(CFLAGS) -o $@  $^ $(LIBDIR) $(LIBRARIES) 

clean:
	rm -f *.o diag

.c.o : ; $(COMPILER) $(CFLAGS) -c $<

#include make.inc

