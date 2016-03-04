#include <stdlib.h>
#include "dspline.h"

void dsp_dscale_4d(dspline *dsp, double *r)
{
	int i, j, k, l;

	double x, y, z, w;
	int tmp;
	double tmptmp;

	int *n = dsp->n;
	double *alpha = dsp->alpha;

	for(l=0;l<n[3];l++)
	{
		w = alpha[3]*2.0;
		if(l == 0)     {w-=alpha[3];}
		if(l == n[3]-1){w-=alpha[3];}
		for(k=0;k<n[2];k++)
		{
			z = alpha[2]*2.0;
			if(k == 0)     {z-=alpha[2];}
			if(k == n[2]-1){z-=alpha[2];}
			for(j=0;j<n[1];j++)
			{
				y = alpha[1]*2.0;
				if(j == 0)     {y-=alpha[1];}
				if(j == n[1]-1){y-=alpha[1];}
				for(i=0;i<n[0];i++)
				{
					tmp = l * n[2] * n[1] * n[0] + k * n[1] * n[0] + j * n[0] + i;
					x = alpha[0]*2.0;
					if(i == 0)     {x-=alpha[0];}
					if(i == n[0]-1){x-=alpha[0];}

					tmptmp = w + z + y + x;
		
					r[tmp] = 1.0 / (w*w + z*z + y*y + x*x + tmptmp*tmptmp) * r[tmp];
				}
			}
		}
	}
}

void dsp_dscale_3d(dspline *dsp, double *r)
{
	int i, j, k;

	double x, y, z;
	int tmp;
	double tmptmp;

	int *n = dsp->n;
	double *alpha = dsp->alpha;

	for(k=0;k<n[2];k++)
	{
		z = alpha[2]*2.0;
		if(k == 0)     {z-=alpha[2];}
		if(k == n[2]-1){z-=alpha[2];}
		for(j=0;j<n[1];j++)
		{
			y = alpha[1]*2.0;
			if(j == 0)     {y-=alpha[1];}
			if(j == n[1]-1){y-=alpha[1];}
			for(i=0;i<n[0];i++)
			{
				tmp = k * n[1] * n[0] + j * n[0] + i;
				x = alpha[0]*2.0;
				if(i == 0)     {x-=alpha[0];}
				if(i == n[0]-1){x-=alpha[0];}

				tmptmp = z + y + x;
	
				r[tmp] = 1.0 / (z*z + y*y + x*x + tmptmp*tmptmp) * r[tmp];
			}
		}
	}
}

void dsp_dscale_2d(dspline *dsp, const double *r, double *z)
{
	int i, j;

	double xs, xe, ys, ye;
	int tmp;
	double tmptmp;

	int *n = dsp->n;
	double *alpha = dsp->alpha;

	for(j=0;j<n[1];j++)
	{
		ys = alpha[1];
		ye = alpha[1];
		if(j == 0)     {ys=0.0;}
		if(j == n[1]-1){ye=0.0;}
		for(i=0;i<n[0];i++)
		{
			tmp = j * n[0] + i;
			xs = alpha[0];
			xe = alpha[0];
			if(i == 0)     {xs=0.0;}
			if(i == n[0]-1){xe=0.0;}

			tmptmp = (ys+ye) + (xs+xe);

			z[tmp] = 1.0 / (ys*ys + ye*ye + xs*xs + xe*xe + tmptmp*tmptmp) * r[tmp];
		}
	}
}

void dsp_stencil_4d(dspline *dsp, const double *x, double *y)
{
	int i, j, k, l;

	double xs, xe, ys, ye, zs, ze, ws, we;
	int tmp;

	int *n = dsp->n;
	double *alpha = dsp->alpha;


	for(l=0;l<dsp->n[3];l++)
	{
		ws = 1.0;
		we = 1.0;
		if(l == 0)     {ws = 0.0;}
		if(l == n[3]-1){we = 0.0;}
		for(k=0;k<n[2];k++)
		{
			zs = 1.0;
			ze = 1.0;
			if(k == 0)     {zs = 0.0;}
			if(k == n[2]-1){ze = 0.0;}
			for(j=0;j<n[1];j++)
			{
				ys = 1.0;
				ye = 1.0;
				if(j == 0)     {ys = 0.0;}
				if(j == n[1]-1){ye = 0.0;}
				for(i=0;i<n[0];i++)
				{
					tmp = l * n[0] * n[1] * n[2] + k * n[0] * n[1] + j * n[0] + i;
					xs = 1.0;
					xe = 1.0;
					if(i == 0)     {xs = 0.0;}
					if(i == n[0]-1){xe = 0.0;}
		
					y[tmp] = ws * alpha[3] * x[tmp - (int)ws * n[0] * n[1] * n[2]]
					       + zs * alpha[2] * x[tmp - (int)zs * n[0] * n[1]]
					       + ys * alpha[1] * x[tmp - (int)ys * n[0]]
					       + xs * alpha[0] * x[tmp - (int)xs]
					       + (-(ws+we) * alpha[3] + -(zs+ze) * alpha[2] + -(ys+ye) * alpha[1] + -(xs+xe) * alpha[0]) * x[tmp]
					       + xe * alpha[0] * x[tmp + (int)xe]
					       + ye * alpha[1] * x[tmp + (int)ye * n[0]]
					       + ze * alpha[2] * x[tmp + (int)ze * n[0] * n[1]]
					       + we * alpha[3] * x[tmp + (int)we * n[0] * n[1] * n[2]];
				}
			}
		}
	}
}

void dsp_stencil_3d(dspline *dsp, const double *x, double *y)
{
	int i, j, k;

	double xs, xe, ys, ye, zs, ze;
	int tmp;

	int *n = dsp->n;
	double *alpha = dsp->alpha;


	for(k=0;k<n[2];k++)
	{
		zs = 1.0;
		ze = 1.0;
		if(k == 0)     {zs = 0.0;}
		if(k == n[2]-1){ze = 0.0;}
		for(j=0;j<n[1];j++)
		{
			ys = 1.0;
			ye = 1.0;
			if(j == 0)     {ys = 0.0;}
			if(j == n[1]-1){ye = 0.0;}
			for(i=0;i<n[0];i++)
			{
				tmp = k * n[0] * n[1] + j * n[0] + i;
				xs = 1.0;
				xe = 1.0;
				if(i == 0)     {xs = 0.0;}
				if(i == n[0]-1){xe = 0.0;}
	
				y[tmp] = zs * alpha[2] * x[tmp - (int)zs * n[0] * n[1]]
				       + ys * alpha[1] * x[tmp - (int)ys * n[0]]
				       + xs * alpha[0] * x[tmp - (int)xs]
				       + (-(zs+ze) * alpha[2] + -(ys+ye) * alpha[1] + -(xs+xe) * alpha[0]) * x[tmp]
				       + xe * alpha[0] * x[tmp + (int)xe]
				       + ye * alpha[1] * x[tmp + (int)ye * n[0]]
				       + ze * alpha[2] * x[tmp + (int)ze * n[0] * n[1]];
			}
		}
	}
}

void dsp_stencil_2d(dspline *dsp, const double *x, double *y)
{
	int i, j;

	double xs, xe, ys, ye;
	int tmp;

	int *n = dsp->n;
	double *alpha = dsp->alpha;


	for(j=0;j<n[1];j++)
	{
		ys = 1.0;
		ye = 1.0;
		if(j == 0)     {ys = 0.0;}
		if(j == n[1]-1){ye = 0.0;}
		for(i=0;i<n[0];i++)
		{
			tmp = j * n[0] + i;
			xs = 1.0;
			xe = 1.0;
			if(i == 0)     {xs = 0.0;}
			if(i == n[0]-1){xe = 0.0;}

			y[tmp] = ys * alpha[1] * x[tmp - (int)ys * n[0]]
			       + xs * alpha[0] * x[tmp - (int)xs]
			       + (-(ys+ye) * alpha[1] + -(xs+xe) * alpha[0]) * x[tmp]
			       + xe * alpha[0] * x[tmp + (int)xe]
			       + ye * alpha[1] * x[tmp + (int)ye * n[0]];
		}
	}
}

void dsp_error_ete(dspline *dsp, const double *x, double *y)
{
	int i;
	int *ytof = dsp->ytof;
	int *used = dsp->used;

	for(i=0;i<dsp->ud;i++)
	{
		y[ytof[used[i]]] = x[ytof[used[i]]];
	}
}

void dsp_error(dspline *dsp, const double *x, double *y)
{
	int i;
	int *ytof = dsp->ytof;
	int *used = dsp->used;

	for(i=0;i<dsp->ud;i++)
	{
		y[i] = x[ytof[used[i]]];
	}
}

void dsp_error_t(dspline *dsp, const double *x, double *y)
{
	int i;
	int *ytof = dsp->ytof;
	int *used = dsp->used;

	for(i=0;i<dsp->ud;i++)
	{
		y[ytof[used[i]]] = x[i];
	}
}

void dsp_gs_2d(dspline *dsp, double *x, const double *b)
{
	int i, j;

	double xs, xe, ys, ye;
	int tmp;
	double tmptmp;

	int *n = dsp->n;
	double *alpha = dsp->alpha;


	for(j=0;j<n[1];j++)
	{
		ys = 1.0;
		ye = 1.0;
		if(j == 0)     {ys = 0.0;}
		if(j == n[1]-1){ye = 0.0;}
		for(i=0;i<n[0];i++)
		{
			tmp = j * n[0] + i;
			xs = 1.0;
			xe = 1.0;
			if(i == 0)     {xs = 0.0;}
			if(i == n[0]-1){xe = 0.0;}

			tmptmp = ys * alpha[1] * x[tmp - (int)ys * n[0]]
			       + xs * alpha[0] * x[tmp - (int)xs]
			       + xe * alpha[0] * x[tmp + (int)xe]
			       + ye * alpha[1] * x[tmp + (int)ye * n[0]];
			x[tmp] = 1.0 / (-(ys+ye) * alpha[1] + -(xs+xe) * alpha[0]) * (b[i] - tmptmp);
		}
	}
}

void dsp_gs_kernel(dspline *dsp, double *x, const double *b)
{
	int i;
	double *y;

	y = malloc(sizeof(double) * (dsp->nn+dsp->ud));

	dsp_stencil_2d(dsp, x, y);
	dsp_error(dsp, x, &y[dsp->nn]);
	dsp_stencil_2d(dsp, y, y);
	dsp_error_t(dsp, &y[dsp->nn], y);
	for(i=0;i<dsp->nn;i++)
	{
		y[i] = b[i] - y[i];
	}
	dsp_dscale_2d(dsp, y, y);

	for(i=0;i<dsp->nn;i++)
	{
		x[i] += y[i];
	}
	free(y);
}
