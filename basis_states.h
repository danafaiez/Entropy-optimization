#ifndef BASIS_STATES_H
#define BASIS_STATES_H

#include <stdio.h>
#include "fock.h"
#include "data_structures.h"
#include "coarse_grain.h"

long long int choose (int N, int Ne);
void  print_binary(int n, int numsites);
void  print_binary_(FILE * file, int n, int numsites);
int iterate_r_basis(int num_sites, int num_particles, int pos[]);
void high_to_low(int num_sites, int num_particles, int * pos, int * reordered);
unsigned long positions_to_binary(int * pos, int num_particles);
ull * enumerate_r_basis(int num_sites, int num_particles);
int binary_to_1d_positions(int * one_d_positions, ull binary);
void one_d_to_2d(int Lx, COORD *c, int * oned_position, int n);
int num_ones(long n);//compute number of ones in binary representation of an integer, e.g. num_ones(33) = 2.
int num_ones_in_range(int start,int end, long n);//compute number of ones in binary representation of an integer, e.g. num_ones(33) = 2,
double L2(double * arr, int length);
HTABLE * create_hash(ull * states);
ull ** calc_regions(PARAMS * pm);//for the moment, the special case of just two separate regions in Hilbert space
ull * init_bases(int num_bath_sites, int num_sites, int num_particles);
double ObsEntropyXE(PARAMS * pm, CG * cg, double * psiE, double * evalues, _Complex double * coef, double t);
double ObsEntropyEX(PARAMS * pm, CG * cg, double * psiEs, double * evalues, _Complex double * psi);
void json_print_density(PARAMS * pm,CG * cg,FILE * f);
_Complex double * psi_thermal(int  numstates,double beta, double * energy, double * evectors, double * e);
_Complex double * coeff(_Complex double * psi, int n, double * h);
double ObsEntropyEX_micro(PARAMS * pm, CG * cg, double * psiEs, double * evalues,double ** evectors);

#endif//BASIS_STATES_H
