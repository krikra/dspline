#include <stdio.h>
#include <stdlib.h>
#include "ippe.h"
#include "ippe_sys.h"
#include "ippe_id.h"
#include "dsp_set.h"

void *f_ippe_init_2d_(int *n1, int *n2)
{
	char lhd[256];

	printf("allocating IPPE\n");

	void *p = malloc(sizeof(ippe));
	ippe *ip = (ippe *)p;
	dspline *dsp;

	printf("setting IPPE parameters\n");

	ip->dim = 2;
	ip->bestv = -1.e+5;
	ip->nextp = 0;
	ip->bestp = 0;
	ip->prebestp = 0;
	ip->iter = 0;
	ip->count = 1;
	ip->terminated = 0;
	ip->term = 5;
	ip->mm = MIN;

	ip->nd_whole = malloc(sizeof(int) * 2);
	ip->nd_whole[0] = *n1;
	ip->nd_whole[1] = *n2;
	ip->dd_whole = (*n1) * (*n2);
	ip->dd_init = 16;

	ip->id = calloc(ip->dd_init, sizeof(int));
	ip->used = calloc(ip->dd_whole, sizeof(int));

	printf("generating initial-design\n");

	sprintf(lhd, "id/lhd/lhd_%d_%d.csv", 16, 1);

	initial_user_2d(ip, lhd);

	printf("setting d-Spline\n");

	ip->dsp = malloc(sizeof(dspline));
	dsp = ip->dsp;

	dsp->itvl = E;

	dspline_setdims(dsp, ip->dim);

	dsp->nd[0] = *n1;
	dsp->nd[1] = *n2;
	dsp->dd = (*n1) * (*n2);

	dsp->n[0] = (*n1) * 3 - 2;
	dsp->n[1] = (*n2) * 3 - 2;
	dsp->nn = dsp->n[0] * dsp->n[1];

	dsp->alpha[0] = 1.e-2;
	dsp->alpha[1] = 1.e-2;

	dspline_setvecs(dsp);

	dspline_set_ytof(dsp);
	
	dspline_set_itr_left(ip->dsp);

	printf("initialize end\n");

	return(p);
}

void f_ippe_choose_(void **ip, int *para)
{
	ippe_choose(*(ippe **)ip, para);
}

void f_ippe_update_(void **ip, int *para, double *val)
{
	ippe_update(*(ippe **)ip, *para, *val);
}

int f_ippe_terminated_(void **ip)
{
	ippe *tmp = *(ippe **)ip;
	return(tmp->terminated);
}

int f_ippe_iteration_(void **ip)
{
	ippe *tmp = *(ippe **)ip;
	return(tmp->iter);
}

void f_ippe_destroy_(void **ip)
{
	int i;
	ippe *tmp = *(ippe **)ip;
	free(tmp->dsp->slv);
	free(tmp->dsp->y);
	free(tmp->dsp->b);
	free(tmp->dsp->f);
	free(tmp->dsp->used);
	free(tmp->dsp->ytof);
	for(i=0;i<tmp->dsp->dim;i++)
	{
		free(tmp->dsp->d_ind[i]);
	}
	free(tmp->dsp->d_ind);
	free(tmp->dsp->alpha);
	free(tmp->dsp->n);
	free(tmp->dsp->nd);
	free(tmp->dsp);

	free(tmp->used);
	free(tmp->id);
	free(tmp);
}
