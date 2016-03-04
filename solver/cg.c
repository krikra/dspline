#include <stdio.h>
#include <stdlib.h>

#include "solver.h"
#include "precon.h"

#include "myblas.h"

#define PASS printf("pass\n");

void cg_kernel(solver *slv)
{
	int i;
	crs *A;
	double *x;
	double *b;
	double *rhist;
	int n;

	int itr;
	double *w;
	double *r,*z,*q,*p,*rr;

	double rho,rho_old,alpha,beta,pq;
	double res,res_init;

	double *old_r;
	double tmp;

	A=slv->A;
	x=slv->x;
	b=slv->b;
	n=slv->n;

	r = calloc(n, sizeof(double));
	z = calloc(n, sizeof(double));
	q = calloc(n, sizeof(double));
	p = calloc(n, sizeof(double));
	//rr = calloc(n, sizeof(double));

	A->spmv(A,x,r);

	xpay(-1.0,b,r,n);

	fill(0.0,p,n);

	rho_old=1.0;

	rho=0.0;

	alpha=0.0;

	beta=0.0;

	norm(&res_init,b,n);

	printf("RESIDUAL 0 %e\n",res_init);

	for(itr=0;itr<slv->mxitr;itr++)
	{
		slv->precon(A,z,r,n);

		dot(&rho,r,z,n);

		beta = rho / rho_old;

		xpay(beta,z,p,n);

		A->spmv(A,p,q);

		dot(&pq,p,q,n);

		alpha = rho / pq;

		axpy(alpha,p,x,n);

		axpy(-alpha,q,r,n);

/*
		A->spmv(A,x,rr);
		xpay(-1.0,b,rr,n);
*/

		norm(&res,r,n);

		res = res / res_init;
		printf("RESIDUAL %d %e\n",itr+1,res);
		if(res<=slv->tol)
		{
			slv->itr=itr+1;
			printf("REACH TO TOLERANCE\n");
			printf("ITERATIONS:%d\n",itr+1);
			break;
		}
		rho_old = rho;
	}

	free(r);
	free(z);
	free(q);
	free(p);
}
