//#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "crs.h"
#include "dspline.h"

#define TRI_RANGE  (2 * range + 1)

void givens(double c, double s, double *x, double *y, int range)
{
	int i;
	double xx,yy;

	for(i=0;i<range;i++)
	{
		xx = x[i]*c + y[i]*s;
		yy =-x[i]*s + y[i]*c;

		x[i] = xx;
		y[i] = yy;
	}
}

/*
void init_mtx(crs *Z)
{
	Z = malloc(sizeof(crs));
	
	Z->val = (double *)calloc((dsp->nn*TRI_RANGE), sizeof(double));
}
*/

void dspline_qr_init(dspline *dsp)
{
	double *work;
	double *Z;
	double r,c,s;
	int i,j,k;

	int bnd[2];

	int coo,tmp,tmptmp;

	int *n;
	int nn;
	int range;
	int dim;

	int hoge;
	int jj;

	n = dsp->n;
	nn = dsp->nn;
	dim = dsp->dim;
	range = nn / n[dsp->dim-1];

	//dsp->Z = alloc_DNS(dsp->nn,TRI_RANGE);

	Z = dsp->Z->val;

	work=malloc(sizeof(double)*(TRI_RANGE));
	//work=malloc(sizeof(double)*(nn));

	for(i=0;i<nn;i++)
	{
		//printf("i:%d\n",i);
		//clean workvector
		for(j=0;j<TRI_RANGE;j++)
		{
			work[j] = 0.0;
		}

		//struct D_i
		tmp = i;
		tmptmp = nn;
		for(j=dim-1;j>=0;j--)
		{
			tmptmp /= n[j];
			coo = tmp / tmptmp;
			tmp = tmp % tmptmp;
			bnd[0] = 0; bnd[1] = 0;
			if(coo-1 >= 0){bnd[0] = 1;}
			if(coo+1 < n[j]){bnd[1] = 1;}
			//printf("bnd:%d %d\n",bnd[0],bnd[1]);
			//if(i - bnd[0] * tmptmp < 0){printf("OOPS \n");}
			//if(i + bnd[1] * tmptmp >= nn){printf("OOPS\n");}
			work[range - bnd[0] * tmptmp] += bnd[0] * dsp->alpha[j];
			work[range + bnd[1] * tmptmp] += bnd[1] * dsp->alpha[j];
			work[range] -= (bnd[0] + bnd[1]) * dsp->alpha[j];
		}

		//givens
		for(j=-range;j<1;j++)
		{
			jj = j + range;
			if(work[jj]==0.0){/*printf("cont\n");*/continue;}
			//printf("%d %d %d %d\n",i,j,jj,TRI_RANGE-jj);
			c = Z[(i+j)*(TRI_RANGE)];
			s = work[jj];
			r = sqrt(c*c + s*s);
			c /= r;
			s /= r;
			//hoge = j + TRI_RANGE >= nn ? nn-j : TRI_RANGE ;
			givens(c,s,&Z[(i+j)*(TRI_RANGE)],&work[jj],TRI_RANGE-jj);
		}
/*
		for(j=0;j<nn;j++)
		{
			printf("%e\n",work[j]);
		}
*/
	}
	free(work);
}

double dspline_qr_incr(dspline *dsp, int ind, double val)
{
	int i,j,k;
	double r,c,s;
	double *work;
	double tmp[2];
	double *Z;

	double res;
	int pos;


	int hoge;

	int *n;
	int nn;
	int range;
	n = dsp->n;
	nn = dsp->nn;
	Z = dsp->Z->val;
	range = nn / n[dsp->dim-1];
	
	//res=dsp->y[ind];
	res = val;
	pos = dsp->ytof[ind];

	work=malloc(sizeof(double)*nn);

	for(i=0;i<nn;i++)
	{
		work[i] = 0.0;
	}
	work[pos] = 1.0;

	for(j=pos;j<nn;j++)
	{
		if(work[j]==0.0){continue;}
		c = Z[j*(TRI_RANGE)];
		s = work[j];
		r = sqrt(c*c + s*s);
		c /= r;
		s /= r;
		hoge = j + TRI_RANGE >= nn ? nn-j : TRI_RANGE ;
		givens(c,s,&Z[j*(TRI_RANGE)],&work[j],hoge);

		tmp[0] = dsp->b[j];
		tmp[1] = res;

		dsp->b[j] = tmp[0]*c + tmp[1]*s;
		res = -tmp[0]*s + tmp[1]*c;
	}
	//res=res;
	free(work);
	return(res);
}

void dspline_backward(dspline *dsp)
{
	int i,j; 

	double *Z,*f,*b;
	int *n;

	int js;
	int nn;
	int range;

	int je;

	double tmp;

	Z = dsp->Z->val;
	f = dsp->f;
	b = dsp->b;
	n = dsp->n;
	nn = dsp->nn;

	range = 2 * nn / n[dsp->dim-1] + 1;

	for(i=nn-1;i>=0;i--)
	{
		tmp = b[i];
		//j = nn-1;
		//while(j>i)
		je = i + range > nn ? nn - i : range ;
		//js = i+1;
		for(j=1;j<je;j++)
		{
			tmp -= Z[i*range+j] * f[i+j];
		}
		f[i] = tmp / Z[i*range];
	}
}
