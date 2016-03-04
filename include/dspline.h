//#include "crs.h"
#include "solver.h"

#ifndef __MY_DSPLINE__
#define __MY_DSPLINE__

typedef enum{
	E,
	X
}interval_strategy;

typedef enum{
	DIR,
	ITR
}dsp_solver;

struct _dspline{
	int dim;
	int *n, *nd;
	int nn, dd, ud;

	int **d_ind;
	int *ytof;

	int *used;
	//int *used_ind;

	interval_strategy itvl;

	int ext;

	double alpha_base;
	double *alpha;

	crs *Z;
	double *f;
	double *b;

	double **x;
	double *y;

	//dsp_Solver s;
	solver *slv;

	void (*solve)(struct _dspline *);
	void (*update)(struct _dspline *, int, double);
};

typedef struct _dspline dspline;

#endif
