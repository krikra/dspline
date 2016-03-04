#ifndef my_MTX
#define my_MTX

struct _crs{
	double *val;
	int *ind;
	int *ptr;
	int nnz;
	int n;

	int symm;
	int diag;

	void (*spmv)(const struct _crs *, const double *, double *);
};

typedef struct _crs crs;

#endif
