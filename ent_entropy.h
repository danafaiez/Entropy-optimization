#ifndef ENT_ENTROPY_H
#define ENT_ENTROPY_H

#include <complex.h>

double ** red_dens_mat(int num_sites, int num_particles, int num_bath_sites, unsigned long * binary_basis, unsigned  long nbasis, double * evector);
//double ** calc_ent_entropy(char * evector_in_name, PARAMS * pm);
double ** calc_ent_entropy(double * evectors, PARAMS * pm);
double calc_ent_entropy_one_ev_(double * evector, PARAMS * pm, int num_bath_sites);
double calc_ent_entropy_one_ev_complex_(_Complex double * evector, PARAMS * pm, int num_bath_sites);
double * calc_ent_entropy_one_ev(double * evector, PARAMS * pm);
double * calc_ent_entropy_n_evs(double ** evectors, PARAMS * pm, int n);
double * calc_ent_entropy_n_evs_superpose(double ** evectors, PARAMS * pm, int n);
double calc_dom_entropy_block(double * state_vector, PARAMS * pm);
double calc_dom_entropy_regions(double * state_vector, PARAMS * pm, ull ** regions);
#endif// ENT_ENTROPY_H
