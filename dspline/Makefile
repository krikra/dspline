CC=gcc

include ../Makefile.in

#INC=../include
#CFLAGS= -O0 -lm -fopenmp -g

OBJ=dsp_abic.o dsp_set.o dsp_qr.o dsp_neq.o dsp_neq_solver.o dsp_stencil.o dsp_ind.o dsp_prob.o dsp_mem.o

dspline.a: ${OBJ}
	ar rcs $@ $^

%.o: %.c
	${CC} ${CFLAGS} -c -o $@ $^ -I${INC}

clean: ;
	rm -f *.a
	rm -f *.o
