#ifndef LOCATE_MULTIVS_H
#define LOCATE_MULTIVS_H
#include "data_structures.h"
#include "basis_states.h"
#include "defs.h"
#include "io.h"
#include "hist.h"

void hist_evector(PARAMS * pm, double * evector, int is_cmplx);
int * top_psi_indices(PARAMS *pm, double * evector, int top_n, int is_cmplx);

#endif//LOCATE_MULTIVS_H
