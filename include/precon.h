#include "solver.h"

#ifndef __MY_PRECON__
#define __MY_PRECON__

void precon_none(const crs *, double *, const double *, const int);
void precon_jacobi_dia(const crs *, double *, const double *, const int);
void precon_sgs_dia(const crs *, double *, const double *, const int);

#endif
