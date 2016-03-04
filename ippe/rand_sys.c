#include <stdlib.h>
#include <stdio.h>
#include <float.h>
#include <math.h>

#include "ippe.h"

#define PASS printf("pass\n");

int random_select(ippe *ip)
{
	int tmp;

	int *used = ip->used;
	int dd = ip->dd_whole;

	do
	{
		tmp = (int)((double)(rand())/(double)(RAND_MAX) * dd);
	}while(used[tmp]);

	return(tmp);
}

void rand_choose(ippe *ip, int *para)
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
		//ip->nextp = ip->id[iter];
		ip->nextp = random_select(ip);
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

void rand_update(ippe *ip, int para, double val)
{
	double tmp;
	if(ip->terminated)
	{
		printf("estimation completed; pass-through\n");
		goto end;
	}
	else
	{
		ip->used[para] = 1;
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
