#include <stdio.h>
#include <stdlib.h>

#include "dsp_stencil.h"
#include "myblas.h"

void dspline_neq_cgnr_stcl(dspline *dsp)
{
	int n;

	int itr;

	double *f = dsp->f;
	double *b = dsp->b;

	double *r, *z, *p, *ap, *q;

	double alpha, beta, rho, rho_old, pq;
	double res, res_init;

	void (*dsp_stencil)(dspline *, const double *, double *);

	switch(dsp->dim)
	{
		case 2:
			dsp_stencil = dsp_stencil_2d; break;
		case 3:
			dsp_stencil = dsp_stencil_3d; break;
		case 4:
			dsp_stencil = dsp_stencil_4d; break;
		default:
			printf("illegal dimension\n");
	}

	r = calloc(dsp->nn, sizeof(double));
	z = calloc(dsp->nn, sizeof(double));
	p = calloc(dsp->nn, sizeof(double));
	ap = calloc(dsp->nn+dsp->ud, sizeof(double));
	q = calloc(dsp->nn, sizeof(double));

	//dsp_stencil(dsp, f, r, ap);
	dsp_stencil(dsp, f, ap);
	dsp_stencil(dsp, ap, r);
	dsp_error_ete(dsp, f, r);

	fill(0.0, ap, n);

	xpay(-1.0, b, r, dsp->nn);

	dsp_gs_2d(dsp, z, r);

	rho_old = 1.0;

	rho = 0.0;

	alpha = 0.0;

	beta = 0.0;

	norm(&res_init, b, dsp->nn);

	printf("RESIDUAL 0 %e\n", res_init);

	for(itr=0;itr<2*dsp->nn;itr++)
	{
		//slv->precon(A,z,r,n);
		copy(r, z, dsp->nn);
		//dsp_gs_2d(dsp, z, r);
		//dsp_dscale_2d(dsp, r, z);

		dot(&rho, r, r, dsp->nn);

		beta = rho / rho_old;

		xpay(beta, z, p, dsp->nn);

		//A->spmv(A,p,q);
		dsp_stencil(dsp, p, ap);
		//dsp_error(dsp, p, &ap[dsp->nn]);
		dsp_stencil(dsp, ap, q);
		//dsp_error_t(dsp, &ap[dsp->nn], q);
		//dot(&pq, ap, ap, dsp->nn+dsp->ud);
		dsp_error_ete(dsp, p, q);
		dot(&pq, p, q, dsp->nn);
		//dsp_gs_2d(dsp, ap, q);

		//ap_dot(&pq, ap, ap, xtnn);

		alpha = rho / pq;

		axpy(alpha, p, f, dsp->nn);

		axpy(-alpha, ap, r, dsp->nn);

		norm(&res, r, dsp->nn);

		res = res / res_init;
		printf("RESIDUAL %d %e\n", itr+1, res);
		if(res <= 1.e-7)
		{
			printf("CONVERGED\n");
			printf("ITERATIONS:%d\n",itr+1);
			break;
		}
		rho_old = rho;
	}
}

void dspline_neq_cgnr(dspline *dsp)
{
	int i;
	int n;
	int itr;

	double *f = dsp->f;
	double *b = dsp->b;

	double *r, *z, *p, *ap, *q, *rr;

	double alpha, beta, rho, rho_old, pq;
	double res, res_init;

	n = dsp->nn + dsp->ud;

	r = calloc(dsp->nn, sizeof(double));
	z = calloc(dsp->nn, sizeof(double));
	p = calloc(dsp->nn, sizeof(double));
	ap = calloc(n, sizeof(double));
	q = calloc(dsp->nn, sizeof(double));
	//rr = calloc(dsp->nn, sizeof(double));

	dsp->Z->spmv(dsp->Z, f, ap);
	//dsp_error(dsp, f, &ap[dsp->nn]);
	dsp->Z->spmv(dsp->Z, ap, r);
	//dsp_error_t(dsp, &ap[dsp->nn], r);
	dsp_error_ete(dsp, f, r);

	fill(0.0, ap, n);

	xpay(-1.0, b, r, dsp->nn);

	//fill(0.0, p, dsp->nn);

	rho_old = 1.0;

	rho = 0.0;

	alpha = 0.0;

	beta = 0.0;

	norm(&res_init, b, dsp->nn);

	printf("RESIDUAL 0 %e\n", res_init);

	for(itr=0;itr<2*dsp->nn;itr++)
	{
		//slv->precon(A,z,r,n);
		copy(r, z, dsp->nn);

		dot(&rho, r, z, dsp->nn);

		beta = rho / rho_old;

		xpay(beta, z, p, dsp->nn);

		//A->spmv(A,p,q);
		dsp->Z->spmv(dsp->Z, p, ap);
		//dsp_error(dsp, p, &ap[dsp->nn]);
		dsp->Z->spmv(dsp->Z, ap, q);
		//dsp_error_t(dsp, &ap[dsp->nn], q);
		dsp_error_ete(dsp, p, q);
		dot(&pq, p, q, dsp->nn);

		alpha = rho / pq;

		axpy(alpha, p, f, dsp->nn);

		axpy(-alpha, q, r, dsp->nn);

/*
		dsp->Z->spmv(dsp->Z, f, ap);
		dsp_error(dsp, f, &ap[dsp->nn]);
		dsp->Z->spmv(dsp->Z, ap, rr);
		dsp_error_t(dsp, &ap[dsp->nn], rr);
		xpay(-1.0, b, rr, dsp->nn);
*/
		
		norm(&res, r, dsp->nn);

		res = res / res_init;
		printf("RESIDUAL %d %e\n", itr+1, res);
		if(res <= 1.e-7)
		{
			printf("REACH TO TOLERANCE\n");
			printf("ITERATIONS:%d\n",itr+1);
			break;
		}
		rho_old = rho;
	}

	free(r);
	free(z);
	free(p);
	free(ap);
	free(q);
}
