#include <stdio.h>
#include <float.h>
#include <math.h>
#include <omp.h>

#include "ippe.h"

#define PASS printf("pass\n");

int difmax(ippe *ip)
{
	int i,j;

	double *f;
	int *ytof;

	int dim;
	int *nd;
	int dd;
	int *used;

	int tmp,coo,foo;
	double dif;

	int bnd[2];

	double maxv = 0.0;
	int maxp = 0;

	dim = ip->dsp->dim;
	nd = ip->dsp->nd;
	dd = ip->dsp->dd;
	used = ip->used;
	f = ip->dsp->f;
	ytof = ip->dsp->ytof;

	for(i=0;i<dd;i++)
	{
		if(used[i]){continue;}
		dif = 0.0;
		tmp = dd;
		foo = i;
		for(j=dim-1;j>=0;j--)
		{
			tmp = tmp / nd[j];
			coo = foo / tmp;
			foo = foo % tmp;
			bnd[0] = 0;bnd[1] = 0;
			if(coo-1 >= 0){bnd[0]=1;}
			if(coo+1 < nd[j]){bnd[1]=1;}
			dif += bnd[0] * f[ytof[i-bnd[0]*tmp]];
			dif += bnd[1] * f[ytof[i+bnd[1]*tmp]];
			dif -= (bnd[0] + bnd[1]) * f[ytof[i]];
			dif = fabs(dif);
		}
		if(maxv < dif){maxv = dif;maxp = i;}
	}
	return(maxp);
}

int best(ippe *ip)
{
	int i;

	double *f;
	int *ytof;

	int dd = ip->dsp->dd;

	double tmp;
	int optp;
	double optv;

	f = ip->dsp->f;
	ytof = ip->dsp->ytof;

	optp = 0;
	optv = -DBL_MAX;
	for(i=0;i<dd;i++)
	{
		tmp = ip->mm * f[ytof[i]];
		if(optv < tmp){optv = tmp;optp = i;}
	}
	return(optp);
}

int terminate(ippe *ip)
{
	if(ip->bestp == ip->prebestp)
	{
		ip->count++;
	}
	else
	{
		ip->count = 1;
	}
	printf("ITER %d count %d\n",ip->iter, ip->count);

	if(ip->count >= ip->term)
	{
		printf("estimation end\n");
		ip->terminated = 1;
	}

	ip->prebestp = ip->bestp;

	return(ip->terminated);
}

void ippe_choose(ippe *ip, int *para)
{
	double t;
	int iter = ip->iter;
	dspline *dsp = ip->dsp;

	if(ip->terminated)
	{
		printf("estimation completed; pass-through\n");
		ip->nextp = ip->bestp;
		goto end;
	}
	if(iter < ip->dd_init)
	{
		ip->nextp = ip->id[iter];
	}
	else
	{
		t = omp_get_wtime();
		dsp->solve(dsp);
		t = omp_get_wtime() - t;
		printf("time-solve %e\n", t);
		ip->bestp = best(ip);

		printf("best %d\n",ip->bestp);

		if(terminate(ip)){ip->nextp = ip->bestp;goto end;}

		if(ip->used[ip->bestp])
		{
			ip->nextp = difmax(ip);
		}
		else
		{
			ip->nextp = ip->bestp;
		}
	}

end:;
printf("choose %d\n",ip->nextp);

	*para = ip->nextp;
}

void ippe_update(ippe *ip, int para, double val)
{
	double t;
	if(ip->terminated)
	{
		printf("estimation completed; pass-through\n");
		goto end;
	}
	else
	{
		t = omp_get_wtime();
		ip->dsp->update(ip->dsp, para, val);
		t = omp_get_wtime() - t;
		printf("time-update %e\n", t);
		//ip->dsp->used_ind[ip->iter] = para;
		ip->used[para] = 1;
		ip->iter++;
	}

end:;
}
