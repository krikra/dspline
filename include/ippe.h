#include "dspline.h"

#ifndef __MY_IPPE__
#define __MY_IPPE__

typedef enum{
	MIN = -1,
	MAX = 1
}MM;

struct _ippe{
	int dim;

	//int *n;
	//int *nd_init;
	int *nd_whole;

	//int nn;
	int dd_init;
	int dd_whole;

	int iter;
	int count;

	int term;
	int terminated;

	int nextp;
	int prebestp;
	int bestp;
	double bestv;

	int *used;
	int *id;

	MM mm;

	dspline *dsp;

	void (*choose)(struct _ippe *, int *);
	void (*update)(struct _ippe *, int, double);
};

typedef struct _ippe ippe;
#endif
