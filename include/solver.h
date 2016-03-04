#include "crs.h"

#ifndef __MY_SOLVER__
#define __MY_SOLVER__

typedef enum{
	J,
	GS,
	CG,
	CR,
	MINRES,
	AMG
}Kernel;

typedef enum{
	pN,
	pJ,
	pSGS,
	pAMG
}Precon;

struct _solver{
	crs *A;
	double *x;
	double *b;

	int n;

	int itr;
	int mxitr;
	double tol;

	void (*kernel)(struct _solver *);
	void (*precon)(const crs *, double *, const double *, const int n);
};

typedef struct _solver solver;

#endif
