#include <stdlib.h>
#include "dspline.h"

void dspline_setdim(dspline *dsp, int dim)
{
	dsp->dim = dim;
	
	dsp->nd = malloc(sizeof(int) * dim);
	dsp->n = malloc(sizeof(int) * dim);
	dsp->alpha = malloc(sizeof(double) * dim);

	dsp->d_ind = malloc(sizeof(int *) * dim);

	dsp->ud = 0;
}

void dspline_setvec(dspline *dsp)
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
