#ifndef MAX_REGION_PROB_H
#define MAX_REGION_PROB_H
#include "data_structures.h"
#include "basis_states.h"
#include <gsl/gsl_multimin.h>
#include "fock.h"


typedef struct {
                _Complex double * coef; double * psiEs; double * evalues;
                PARAMS * pm;
                CG * cg;
                ull ** regions;

               }PARAM_MAX;
double my_prob_func(const gsl_vector *x, void *params);
ull ** calc_regions_x(PARAMS * pm);
_Complex double * psi_phi_x(const gsl_vector *x, _Complex double * coef, double * psiEs);
int regional_prob_max(PARAMS * pm, CG * cg, _Complex double *coef  ,double * psiEs, double * evalues, ull ** reg);
double L(_Complex double * arr, int n);
//double * den(PARAMS * pm, _Complex double * psi, ull *  binary_basis);
//ull bit_val(ull num,int pos);
//ull bit_val(1&((ull num) >> (int pos)));
#endif//MAX_REGION_PROB
