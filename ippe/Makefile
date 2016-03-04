CC=gcc

include ../Makefile.in

#INC=../include
#CFLAGS= -O3 -lm -fopenmp -g

#SRC=ippe_sys.c 

OBJ=ippe_sys.o ippe_id.o ippe_mem.o bnb_sys.o rand_sys.o

ippe.a: ${OBJ}
	ar rcs $@ $^

%.o: %.c
	${CC} ${CFLAGS} -c -o $@ $^ -I${INC}

clean: ;
	rm -f *.a
	rm -f *.o
