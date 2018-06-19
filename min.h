#ifndef MIN_S_EX_H
#define MIN_S_EX_H
#include <gsl/gsl_multimin.h>
#include "energy_graining.h"


typedef struct {
                _Complex double * coef; double * h; double * psiEs; double * evalues; _Complex double *psi_e_b;    
                PARAMS * pm; EG * eg;
                CG * cg;//take out h later
               }PARAM_MIN;
_Complex double * psi_phi(const gsl_vector *x, _Complex double * coef, double * psiEs);
double my_f_ex(const gsl_vector * x, void * params);
double my_f_ent(const gsl_vector * x, void * params);
double my_f_fac(const gsl_vector * x, void * params);
double Entropy_min(PARAMS * pm, CG * cg, _Complex double *coef, double * psiEs, double * evalues, EG * eg);

double ** calc_PsEX(PARAMS * pm,CG * cg, double * evalues,double * psiEs, _Complex double * psi);
double ** coarse_density(PARAMS * pm, CG *cg, double ** Ps);
double * den(PARAMS * pm, _Complex double * psi, ull *  binary_basis);
#endif//MIN_S_EX_H
