#ifndef THERMO_H
#define THERMO_H

typedef struct {
   double * zarr, * earr;
   double dbeta, betamax;
} THERMO;

THERMO * calc_thermo(double * energy, int N, double dbeta);
int calc_beta_index(THERMO * thermo, double Ein);
double thermal_entropy(THERMO * therm, double Ein, int beta_index);
double S_of_beta(double beta, double * energy, int N);
void Z_of_beta(double beta, double * energy, int N, double * e_out, double * z_out);
double E_of_beta(double beta, double * energy, int N);
double E_of_beta(double beta, double * energy, int N);

#endif
