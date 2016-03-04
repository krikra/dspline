#include <stdlib.h>
#include "ippe.h"
#include "dsp_mem.h"

void ippe_setdim(ippe *ip, int dim)
{
	ip->dim = dim;

	ip->bestv = -1.e+7;
	ip->nextp = 0;
	ip->bestp = 0;
	ip->prebestp = 0;
	ip->iter = 0;
	ip->count = 1;
	ip->terminated = 0;

	ip->nd_whole = malloc(sizeof(int) * dim);
	ip->dsp = malloc(sizeof(dspline));
	dspline_setdim(ip->dsp, ip->dim);
	ip->dsp->itvl = E;
	ip->dsp->ext = 0;
	//ip->n = malloc(sizeof(int) * dim);
}

void ippe_setvec(ippe *ip)
{
	int i;
	ip->id = malloc(sizeof(int) * ip->dd_init);
	ip->used = calloc(ip->dd_whole, sizeof(int));

	ip->dsp->nn = 1;
	for(i=0;i<ip->dim;i++)
	{
		ip->dsp->nd[i] = ip->nd_whole[i];
		ip->dsp->n[i] = ip->nd_whole[i] * 3 - 2;
		ip->dsp->alpha[i] = 1.e-2;
		ip->dsp->nn *= ip->dsp->n[i];
	}
	ip->dsp->dd = ip->dd_whole;
	//ip->dsp->nn = ;

	dspline_setvec(ip->dsp);

	dspline_set_ytof(ip->dsp);
	//dspline_init_mtx(ip->dsp);
}

void ippe_destroy(ippe *ip)
{
	dspline_destroy(ip->dsp);
	free(ip->dsp);
	free(ip->id);
	free(ip->used);
	free(ip->nd_whole);
}
