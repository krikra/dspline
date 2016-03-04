#include <math.h>

void fill(double a,double *x,int n)
{
	int i;
	for(i=0;i<n;i++)
	{
		x[i] = a;
	}
}

void copy(double *x,double *y,int n)
{
	int i;
	for(i=0;i<n;i++)
	{
		y[i] = x[i];
	}
}

void xpay(const double a,const double *x,double *y,const int n)
{
	int i;
	for(i=0;i<n;i++)
	{
		y[i] = x[i] + a * y[i];
	}
}

void axpy(const double a,const double *x,double *y,const int n)
{
	int i;
	for(i=0;i<n;i++)
	{
		y[i] = a * x[i] + y[i];
	}
}

void dot(double *a,const double *x,const double *y,int n)
{
	int i;
	double hoge=0.0;
	for(i=0;i<n;i++)
	{
		hoge = hoge + x[i] * y[i];
	}
	*a = hoge;
}

void norm(double *a,const double *x,int n)
{
	double hoge=0.0;
	dot(&hoge,x,x,n);

	*a = sqrt(hoge);
}
