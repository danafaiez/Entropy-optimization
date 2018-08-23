#ifndef UNITARY_H
#define UNITARY_H
#include "data_structures.h" 
//#define UNITARY_TEST


typedef struct {
                 int N,M;
                 complx ** E, ** newE;
                 double ** J, *a;
                 int numstates, size_of_box;
                 PARAMS * pm;
                 double * psiEs;
                 double * energy;
                 CG * cg;
                 EG * eg;
} PSI_STATE;

typedef struct {
   int i,j;
   PSI_STATE * psi_state;
} PARAM_UNITARY;
//complx ** makeEsN(PARAMS * pm, double * evectors);
ull ** calc_regionsEs(PARAMS * pm);
double unitary_min(CG * cg, PARAMS * pm,  complx *coef, double * psiEs, double *
energy, EG * eg);
complx * make_coeffs(int n,int real_or_complx);
double unitary_test(int N, int M);
complx ** makeEN(PARAMS * pm, double * evectors);
#endif// UNITARY_H
