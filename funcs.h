#ifndef FUNCS_H
#define FUNCS_H
#define _GNU_SOURCE

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>


#include <complex.h>
#include "clapack.h"
#include "clapack.h"
#include "fock.h"
#include "data_structures.h"

long long int choose(int, int);
//void tVmodel_2d(int Ne,int Lx,int Ly,long long int Nbasis,double t,double tprime,double V);
//void tVmodel_2d(PARAMS * pm);
ull a(ull state, int site, int N, int * sign);
ull adag(ull state, int site, int N, int * sign);
long long int binTodec(int *a, int n);
double * diag(double * h, int  numstates, PARAMS * pm);
_Complex double * convert1to2(PARAMS * pm, _Complex double * psi1,int Linit, ull * basis1, int Lfinal, ull * basis2);
_Complex double * psit(_Complex double * coef, double * h, double * evalues, double t);
double * read_energies(char * energy_in_name, PARAMS * pm);
void write_efiles(PARAMS * pm, double * h, double * ener);
FILE * write_energies(FILE * energy_out, double * energies, PARAMS * pm, char * energy_out_name, char * dirname);
_Complex double * convert_double2complex(double * in, int N);
#endif//FUNCS_H
