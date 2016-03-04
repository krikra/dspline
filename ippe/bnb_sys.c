#include <stdio.h>
#include <float.h>
#include <math.h>

#include "ippe.h"

#define PASS printf("pass\n");

void bnb_choose(ippe *ip, int *para)
{
	int iter = ip->iter;
	dspline *dsp = ip->dsp;

	if(ip->terminated)
	{
		printf("estimation completed; pass-through\n");
		//ip->nextp = ip->bestp;
		goto end;
	}
	if(iter < ip->dd_init)
	{
		ip->nextp = ip->id[iter];
	}
	else
	{
		ip->terminated = 1;
		ip->nextp = ip->bestp;
	}

end:;
printf("choose %d\n",ip->nextp);

	*para = ip->nextp;
}

void bnb_update(ippe *ip, int para, double val)
{
	double tmp;
	if(ip->terminated)
	{
		printf("estimation completed; pass-through\n");
		goto end;
	}
	else
	{
		tmp = ip->mm * val;
		if(ip->bestv < tmp)
		{
			ip->bestp = para;
			ip->bestv = tmp;
		}
		ip->iter++;
		//ip->hoge = val;
	}

end:;
}
