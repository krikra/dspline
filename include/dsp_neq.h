#include "dspline.h"

#ifndef __MY_DSP_NEQ__
#define __MY_DSP_NEQ__

void dspline_neq_incr_left(dspline *, int, double);
void dspline_neq_incr_right(dspline *, int, double);
void dspline_neq_init_left(dspline *);
void dspline_neq_init_right(dspline *);
void dspline_neq_solve(dspline *);

#endif
