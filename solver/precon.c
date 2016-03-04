#include "solver.h"

void precon_none(const crs *A, double *z,const double *r,const int n)
{
	int i;
	for(i=0;i<n;i++)
	{
		z[i] = r[i];
	}
}

void precon_jacobi_dia(const crs *A,double *z, const double *r,const int n)
{
	int i;
	for(i=0;i<n;i++)
	{
		z[i]=r[i]/A->val[A->diag*A->n+i];
	}
}

void precon_sgs_dia(const crs *A,double *z, const double *r,const int n)
{
	int i,j,jj;

	double tmp;

	for(i=A->n-1;i>=0;i--)
	{
		tmp = r[i];
		for(j=0;j<A->diag;j++)
		{
			jj=A->ind[j];
			if(i-jj>=n){continue;}
			tmp -= z[i-jj] * A->val[j*A->n+i-jj];
		}
		z[i] = tmp / A->val[A->diag*A->n+i];
	}

	for(i=0;i<A->n;i++)
	{
		z[i] *= A->val[A->diag*A->n+i];
	}

	for(i=0;i<A->n;i++)
	{
		tmp = z[i];
		for(j=0;j<A->diag;j++)
		{
			jj = A->ind[j];
			if(i+jj<0){continue;}
			tmp -= z[i+jj] * A->val[j*A->n+i];
		}
		z[i] = tmp / A->val[A->diag*A->n+i];
	}
}
