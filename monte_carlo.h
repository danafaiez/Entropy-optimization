#ifndef MONTE_CARLO_H
#define MONTE_CARLO_H
#include <gsl/gsl_multimin.h>
#include "energy_graining.h"

double monte_func(PARAMS * pm, CG * cg, _Complex double *coef, double * psiEs, double * evalues, EG * eg, _Complex double * psi_init);

#endif
