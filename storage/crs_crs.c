#include <stdlib.h>
#include <stdio.h>

#include "crs.h"

crs *alloc_crs(int nnz, int n, int symm)
{
	int i;
	crs *A;
	a=malloc(sizeof(crs));
	A->nnz=nnz;
	A->n=n;
	A->val=malloc(sizeof(double)*nnz);
	A->ind=malloc(sizeof(int)*nnz);
	A->ptr=malloc(sizeof(int)*(n+1));
	for(i=0;i<nnz;i++)
	{
		A->val[i]=0.0;
		A->ind[i]=0;
	}
	for(i=0;i<n+1;i++)
	{
		A->ptr[i]=0;
	}

	A->spmv = symm?symv_crs:spmv_crs;

	return(A);
}

void free_crs(crs *A)
{
	free(A->val);
	free(A->ind);
	free(A->ptr);
	free(A);
}

void spmv_crs(const crs *A,const double *x,double *y,const int n)
{
	int i,j;
	for(i=0;i<n;i++)
	{
		y[i]=0.0;
		for(j=A->ptr[i];j<A->ptr[i+1];j++)
		{
			y[i]+=A->val[j]*x[A->ind[j]];
		}
	}
}

void spmtv_crs(const crs *A,const double *x,double *y,const int n)
{
	int i,j;
	for(i=0;i<n;i++)
	{
		y[i]=0.0;
	}
	for(i=0;i<n;i++)
	{
		for(j=A->ptr[i];j<A->ptr[i+1];j++)
		{
			y[A->ind[j]]+=A->val[j]*x[i];
		}
	}
}

void symv_crs(const crs *A,const double *x,double *y,const int n)
{
	int i,j;
	for(i=0;i<n;i++)
	{
		y[i]=0.0;
		for(j=A->ptr[i];j<A->ptr[i+1];j++)
		{
			y[i]+=A->val[j]*x[A->ind[j]];

			y[i + A->ind[j]]+=A->val[j]*x[j];
		}
	}
}

