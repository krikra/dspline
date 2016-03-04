#include <stdio.h>
#include <stdlib.h>

#include "precon.h"
#include "solver.h"

#include "cg.h"

void solver_set_kernel(solver *slv, Kernel kind)
{
	switch(kind)
	{
		case CG:
			slv->kernel = cg_kernel;
			break;
	}
}

void solver_set_precon(solver *slv, Precon kind)
{
	switch(kind)
	{
		case pN:
			slv->precon = precon_none;
			break;

		case pJ:
			slv->precon = precon_jacobi_dia;
			break;

		case pSGS:
			slv->precon = precon_sgs_dia;
			break;
	}
}
