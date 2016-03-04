#include <stdio.h>
#include <math.h>

#include "dspline.h"
#include "dsp_qr.h"

#define TRI_RANGE (2 * range + 1)

static const double golden = (1.0+sqrt(5.0))/2.0;

double get_abic(dspline *dsp,double alpha)
{
	int i;

	double fst = 0.0;
	double trd = 0.0;

	double tmp;

	int nn,dd;

	int range;

	double *Z = dsp->Z->val;

	nn = dsp->nn;
	dd = dsp->dd;
	range = dsp->nn / dsp->n[dsp->dim-1];

	dsp->alpha[0] = alpha;

	dspline_qr_init(dsp);

	for(i=0;i<dd;i++)
	{
		tmp = dspline_qr_incl(dsp,i,dsp->y[i]);
		trd += tmp * tmp;
	}

	for(i=0;i<nn;i++)
	{
		fst+=log(Z[i*(TRI_RANGE)]);
	}
	fst *= 2.0;

	return(fst - 2 * (nn - 2) * log(dsp->alpha[0]) + dd - 2.0 * log(trd));
}

// 4 points -> next pointt & 3 points

int abic_golden(double *key_ext,double *key_int,double *val_int)
{
	int flag;
	//const double golden=(1.0+sqrt(5.0))/2.0;
	flag = val_int[0]>val_int[1]?0:1;

	key_ext[flag] = key_int[flag];

	key_int[flag] = key_int[!flag];
	val_int[flag] = val_int[!flag];

	key_int[!flag] = (1.0*flag+golden*!flag)/(1.0+golden)*fabs(key_ext[1]-key_ext[0])+key_ext[0];

	return(!flag);
}

void dspline_abic_search(dspline *dsp)
{
	int flag;
	//const double golden=(1.0+sqrt(5.0))/2.0;
	double itv=1.0;

	double key_ext[2];
	double key_int[2];
	double val_int[2];

	key_ext[0] = -5.0;
	key_ext[1] = 5.0;

	key_int[0] = 1.0/(1.0+golden)*fabs(key_ext[1]-key_ext[0])+key_ext[0];
	key_int[1] = golden/(1.0+golden)*fabs(key_ext[1]-key_ext[0])+key_ext[0];

	val_int[0] = get_abic(dsp,key_int[0]);
	val_int[1] = get_abic(dsp,key_int[1]);

	while(itv>0.01)
	{
		flag = abic_golden(key_ext, key_int, val_int);
		val_int[flag] = get_abic(dsp, pow(10.0, key_int[flag]));
		itv = fabs((key_int[1]-key_int[0])/key_int[0]);
		printf("alpha=%e ABIC=%e res=%e\n", pow(10.0,key_int[flag]), val_int[flag], itv);
	}
}
