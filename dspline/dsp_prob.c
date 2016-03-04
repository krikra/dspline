//#include <stdio.h>
#include <stdlib.h>
//#include <math.h>

#include "crs_dia.h"

void generate_newD_DIA(crs *D,const int dim,const int *n,const double *alpha,int symm)
{
	int i,j,k,l;
	int nn=1.0;
	for(i=0;i<dim;i++)
	{
		nn*=n[i];
	}

	int tmp;
	int bnd;
	int now;
	int ii;
	int coo;

	tmp=1;
	D->ind[dim]=0;
	D->diag = dim;

	for(i=0;i<dim;i++)
	{
		D->ind[dim-i-1]=-tmp;
		if(!symm){D->ind[dim+i+1]=tmp;}
		tmp*=n[i];
	}

	k=0;
	tmp = nn;
	for(i=0;i<D->diag;i++)
	{
		now = n[D->diag-i-1];
		ii = D->ind[i];
		tmp /= now;
		for(j=-ii;j<D->n;j++)
		{
			coo = j / tmp;
			coo = coo % now;
			bnd = 0;
			if(coo-1>=0){bnd=1;}
			D->val[i*D->n+j]=bnd*alpha[i];
			if(!symm){D->val[(D->nnz-i-1)*D->n+j+ii]=bnd*alpha[i];}
			D->val[D->diag*D->n+j]-=bnd*alpha[i];
			D->val[D->diag*D->n+j+ii]-=bnd*alpha[i];
		}
	}
}

void setind(int *ind,const int dim,const int *n,int interval,const int center)
{
	int i;
	int offset[2];
	int tmp,subcenter;

	offset[1] = dim * dim + dim;
	offset[0] = dim;

	interval /= n[dim-1];

	subcenter = center - offset[1];

	ind[subcenter] = -2 * interval;

	tmp = interval;
	for(i=dim-1;i>0;i--)
	{
		tmp /= n[i-1];
		ind[subcenter+offset[0]-i] = -interval - tmp;
		ind[subcenter+offset[0]+i] = -interval + tmp;
	}
	ind[subcenter+offset[0]] = -interval;

	if(dim != 1){setind(ind,dim-1,n,interval,center);}

	ind[center] = 0;
}

void generate_DtD2_DIA(crs *DtD, const int dim, const int *n, const double *alpha, const int mpi_offset, const int mpi_offset_next, int symm)
{
	int i,j,k;
	crs *D;
	int nn=1;
	for(i=0;i<dim;i++)
	{
		nn*=n[i];
	}
	D=alloc_dia(2*dim+1, nn, 0);
	generate_newD_DIA(D, dim, n, alpha, 0);
	int nnd=D->nnz;

	DtD->diag = dim * dim + dim;

	setind(DtD->ind,dim,n,nn,DtD->diag);

	if(symm == 0)
	{
		for(i=0;i<DtD->nnz/2;i++)
		{
			DtD->ind[DtD->nnz-i-1] = -DtD->ind[i];
		}
	}

	int tmp,tmptmp,tmptmptmp,tmptmptmptmp;
	int offset,unoffset,joffset;
	int is,js,ks;
	int ie,je,ke;
	int udx,ldx,unudx,unldx;
	int where;

	for(j=0;j<nnd;j++)
	{
		tmp=nn-abs(D->ind[j]);
		if(D->ind[j]>=0)
		{
			udx=0;
			ldx=-D->ind[j];
			unudx=D->ind[j];
			unldx=0;
			joffset=0;
		}else
		{
			udx=-D->ind[j];
			ldx=0;
			unudx=0;
			unldx=-D->ind[j];
			joffset=-D->ind[j];
		}
		for(i=0;i<nnd;i++)
		{
			tmptmp=nn-abs(D->ind[i]);
			tmptmptmp=D->ind[i]+D->ind[j];
			if(symm && (tmptmptmp > 0)){continue;}
			for(k=0;k<DtD->nnz;k++)
			{
				if(tmptmptmp==DtD->ind[k]){where=k;break;}
			}
			ks=D->ind[i]>udx?D->ind[i]-udx:0;
			ke=D->ind[i]<ldx?tmp-(ldx-D->ind[i]):tmp;
			offset=D->ind[i]>=udx?-ks:udx-D->ind[i];
			unoffset=DtD->ind[where]>=unudx?-ks:unudx-DtD->ind[where];
			for(k=ks;k<ke;k++)
			{
				DtD->val[where*nn+k+unoffset]+=D->val[i*nn+k+offset]*D->val[j*nn+k+joffset];
			}
		}
	}
	free_dia(D);
}
