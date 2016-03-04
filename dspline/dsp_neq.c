#include <stdio.h>
#include <stdlib.h>

#include "solver.h"
#include "crs_dia.h"
#include "dsp_prob.h"
#include "dspline.h"

void dspline_neq_incr_left(dspline *dsp, int ind, double val)
{
	int i;
	int tmp;
	crs *Z;

	Z = dsp->Z;

	tmp = dsp->ytof[ind];

	Z->val[Z->diag*Z->n+tmp] += 1.0;

	dsp->b[tmp] = val;
	dsp->ud++;
}

void dspline_neq_incr_right(dspline *dsp, int ind, double val)
{
	dsp->used[dsp->ud] = ind;
	dsp->b[dsp->ytof[ind]] = val;
	dsp->ud++;
}

void dspline_neq_init_left(dspline *dsp)
{
	int tmp;
	int symm = 1;
	tmp = dsp->dim * dsp->dim + dsp->dim;

	dsp->Z = alloc_dia(tmp+1, dsp->nn, symm);
	generate_DtD2_DIA(dsp->Z, dsp->dim, dsp->n, dsp->alpha, 0, dsp->nn, symm);
}

void dspline_neq_init_right(dspline *dsp)
{
	dsp->Z = alloc_dia(dsp->dim+1, dsp->nn, 1);
	generate_newD_DIA(dsp->Z, dsp->dim, dsp->n, dsp->alpha, 1);
}

void dspline_neq_solve(dspline *dsp)
{
	solver *slv = dsp->slv;

	slv->b = dsp->b;
	slv->x = dsp->f;
	slv->A = dsp->Z;

	slv->kernel(slv);

	//dsp->Z = slv->A;
	dsp->f = slv->x;
	//dsp->b = slv->b;
}

void dspline_neq_destroy(dspline *dsp)
{
	free_dia(dsp->Z);
}
