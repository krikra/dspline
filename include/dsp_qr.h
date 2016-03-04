#include "dspline.h"

#ifndef __MY_DSP_QR__
#define __MY_DSP_QR__

void dspline_qr_init(dspline *);
double dspline_qr_incr(dspline *, int, double);
void dspline_backward(dspline *);

#endif
