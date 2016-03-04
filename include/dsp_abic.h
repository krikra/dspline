#include "dspline.h"
#include "dsp_qr.h"

#ifndef __MY_DSP_ABIC__
#define __MY_DSP_ABIC__

double get_abic(dspline *dsp);
int abic_golden(dspline *dsp);
void dspline_abic_search(dspline *dsp);

#endif
