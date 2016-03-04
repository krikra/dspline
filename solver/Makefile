CC=gcc

include ../Makefile.in
#INC=../include
#CFLAGS= -O3 -lm -fopenmp -g

solver.a: solver.o cg.o precon.o myblas.o
	ar rcs $@ $^

%.o: %.c
	${CC} ${CFLAGS} -c -o $@ $^ -I${INC}

clean: ;
	rm -f *.a
	rm -f *.o
