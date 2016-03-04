#include <stdlib.h>

#include "dsp_qr.h"
#include "dsp_neq.h"
#include "dsp_neq_solver.h"

/*
void dspline_setdims(dspline *dsp, int dim)
{
	dsp->dim = dim;
	
	dsp->nd = malloc(sizeof(int) * dim);
	dsp->n = malloc(sizeof(int) * dim);
	dsp->alpha = malloc(sizeof(double) * dim);

	dsp->d_ind = malloc(sizeof(int *) * dim);

	dsp->ud = 0;
}

void dspline_setvecs(dspline *dsp)
{
	int i;
	
	for(i=0;i<dsp->dim;i++)
	{
		dsp->d_ind[i] = malloc(sizeof(int) * dsp->nd[i]);
	}
	dsp->ytof = calloc(dsp->dd, sizeof(int));
	dsp->used = calloc(dsp->dd, sizeof(int));

	dsp->f = calloc(dsp->nn, sizeof(double));
	dsp->b = calloc(dsp->nn, sizeof(double));
	dsp->y = calloc(dsp->dd, sizeof(double));
}
*/

void dspline_set_dir(dspline *dsp)
{
	dsp->update = dspline_qr_incr;
	dsp->solve = dspline_backward;

	dspline_qr_init(dsp);
}

void dspline_set_itr_left(dspline *dsp)
{
	dsp->update = dspline_neq_incr_left;
	dsp->solve = dspline_neq_solve;

	dspline_neq_init_left(dsp);

	dsp->slv = malloc(sizeof(solver));
	solver_set_kernel(dsp->slv, CG);
	solver_set_precon(dsp->slv, pSGS);
	dsp->slv->n = dsp->nn;
	dsp->slv->mxitr = 2 * dsp->nn;
	dsp->slv->tol = 1.e-7;
}

void dspline_set_itr_right(dspline *dsp)
{
	dsp->update = dspline_neq_incr_right;
	dsp->solve = dspline_neq_cgnr;

	dspline_neq_init_right(dsp);
}

void dspline_set_itr_right_stcl(dspline *dsp)
{
	dsp->update = dspline_neq_incr_right;
	dsp->solve = dspline_neq_cgnr_stcl;
}

/*
void dspline_destroy(dspline *dsp)
{
	int i;

	free(dsp->y);
	free(dsp->b);
	free(dsp->f);

	free(dsp->used);
	free(dsp->ytof);

	for(i=0;i<dsp->dim;i++)
	{
		free(dsp->d_ind[i]);
	}
	free(dsp->d_ind);

	free(dsp->alpha);
	free(dsp->n);
	free(dsp->nd);
}
*/
