#include <stdio.h>
#include <math.h>
#include "defs.h"
#include "thermo.h"
//#include "htable.h"
//#include "funcs.h"
//#include "fock.h"
//#include "basis_states.h"
//#include "coarse_grain.h"

void Z_of_beta(double beta, double * energy, int N, double * e_out, double * z_out)
{
      int i;
      double z=0;
      double e, e_num = 0;
      for(i=0;i < N;i++)
      {
         z += exp(-beta*energy[i]);
         e_num  += energy[i]*exp(-beta*energy[i]);
      }
      *e_out = e_num/z;
      *z_out = z;
}


double E_of_beta(double beta, double * energy, int N)
{
   double e, z;
   Z_of_beta(beta, energy, N,  &e, &z);
   return e;
}

double S_of_beta(double beta, double * energy, int N)
{
   double e, z;
   Z_of_beta(beta, energy, N,  &e, &z);
   return beta*e+log(z);
}

THERMO * calc_thermo(double * energy, int N, double dbeta)
{
   THERMO * thermo;
   double * zarr=0;
   double * earr=0;
   double beta;
   double diffE = energy[N-1]-energy[0];
   double betamax = 50.0/diffE;
   for(beta=0; beta < betamax; beta += dbeta)
   {
      double e_out, z_out;
      Z_of_beta(beta, energy, N, &e_out, &z_out);
      appendarr_(zarr,z_out);
      appendarr_(earr,e_out);
   }
   newarr_(thermo,1);
   thermo->zarr = zarr;
   thermo->earr = earr;
   thermo->betamax = betamax;
   thermo->dbeta = dbeta;
   return thermo;
}

int calc_beta_index(THERMO * thermo, double Ein)
{
   double * earr = thermo->earr;
   int n = size_(earr);
   int i_low = 0;
   if (Ein < earr[n-1])
   {
      printf ("Ein is too small\n");
      return -1.0;
   }
   else if (Ein > earr[0])
   {
      printf ("Ein is too large\n");
      return (int) 1e9;
   }
   else
   {
      int imiddle=-1;
      int i_high = size_(earr);

      while (i_high - i_low > 1)
      {
         imiddle = (i_low+i_high)/2;
         if (Ein > earr[imiddle])
         {
            i_high = imiddle;
         }
         else 
         {
            i_low = imiddle;
         }
      }
   }
   return i_low;
}

double thermal_entropy(THERMO * therm, double Ein, int beta_index)
{
   //int beta_index = calc_beta_index(therm, Ein);
   if (beta_index < 9e9 && beta_index >= 0)
      return -1;
   double beta = beta_index*therm->dbeta;
   double T = 1.0/beta;
   double F = -T*log(therm->zarr[beta_index]);
   double S = (Ein-F)/T;
   return S;
}


