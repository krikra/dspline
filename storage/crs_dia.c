#include <stdlib.h>
#include <stdio.h>

#include "crs.h"

void spmv_dia(const crs *A,const double *x,double *y)
{
	int i,j;
	int tmp,offset,uoffset;

	int n = A->n;

	for(i=0;i<n;i++)
	{
		y[i] = 0.0;
	}

	for(i=0;i<A->nnz;i++)
	{
		tmp=A->ind[i];
		if(tmp>=0)
		{
			offset=0;
			uoffset=tmp;
		}
		else
		{
			offset=-tmp;
			uoffset=0;
		}
		tmp=abs(tmp);
		for(j=0;j<A->n-tmp;j++)
		{
			y[j+offset] += A->val[i*A->n+j+offset] * x[j+uoffset];
		}
	}
}

void symv_dia(const crs *A,const double *x,double *y)
{
	int i,j;
	int tmp,offset,uoffset;

	int n = A->n;

	for(i=0;i<n;i++)
	{
		y[i] = 0.0;
	}

	for(i=0;i<A->nnz-1;i++)
	{
		tmp = -A->ind[i];

		for(j=0;j<A->n-tmp;j++)
		{
			y[j + tmp] += A->val[i*A->n + j + tmp] * x[j];

			y[j] += A->val[i*A->n + j + tmp] * x[j + tmp];
		}
	}

	for(j=0;j<A->n;j++)
	{
		y[j] += A->val[(A->nnz-1)*A->n + j] * x[j];
	}
}

crs *alloc_dia(int nnd, int n, int symm)
{
	int i;
	crs *A;
	A=malloc(sizeof(crs));
	A->nnz=nnd;
	A->n=n;
	A->val=malloc(sizeof(double)*nnd*n);
	A->ind=malloc(sizeof(int)*nnd);
	for(i=0;i<nnd*n;i++)
	{
		A->val[i]=0.0;
	}
	for(i=0;i<nnd;i++)
	{
		A->ind[i]=0;
	}

	A->symm = symm;

	A->spmv = symm?symv_dia:spmv_dia;
	//printf("%p %d\n",&A->nnz,A->nnz);fflush(stdout);
	return(A);
}

void free_dia(crs *A)
{
	free(A->val);
	free(A->ind);
	free(A);
}

