#include "dsp_qr.h"
#include "dsp_neq.h"

void dspline_set_dir(dspline *dsp)
{
	dsp->update = dspline_qr_incl;
	dsp->solve = dspline_backward;

	dsp_qr_init(dsp);
}

void dspline_set_itr(dspline *dsp)
{
	dsp->update = dspline_neq_incl;
	dsp->solve = dspline_neq_solve;

	dspline_neq_init(dsp);
}
