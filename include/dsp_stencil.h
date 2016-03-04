#include "dspline.h"

#ifndef __DSP_STENCIL__
#define __DSP_STENCIL__
void dsp_stencil_4d_bb(dspline *, double *, double *);
void dsp_stencil_3d_bb(dspline *, double *, double *);
void dsp_stencil_2d_bb(dspline *, double *, double *);
void dsp_stencil_4d(dspline *, const double *, double *);
void dsp_stencil_3d(dspline *, const double *, double *);
void dsp_stencil_2d(dspline *, const double *, double *);
void dsp_error_ete(dspline *, const double *, double *);
void dsp_error(dspline *, const double *, double *);
void dsp_error_t(dspline *, const double *, double *);
#endif
