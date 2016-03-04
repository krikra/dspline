CC=gcc

include ../Makefile.in
#INC=../include
#CFLAGS= -O3 -fopenmp -g

OBJ= crs_dia.o

storage.a: ${OBJ}
	ar rcs $@ $^

%.o: %.c
	${CC} ${CFLAGS} -c -o $@ $^ -I${INC}

clean: ;
	rm -f *.a
	rm -f *.o
