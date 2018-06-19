 /*this program maximizes the probability that the wavefunction is localized within x_init and x_fin*/ 
#include "locate_multivs.h"
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_vector.h>
#include "max_region_prob.h"
#include "min.h"
#include "ent_entropy.h"

/*unitary evolution of the wavefunction*/
_Complex double * psi_phi_x(const gsl_vector *x, _Complex double * coef, double * psiEs)
 {
    int i;
    int n = size_(coef);
    _Complex double * psi;
    newarr_(psi,n);
    double phi;
    for(i=0;i<n;i++)
    {
       int base = n*i;
       int j;
       phi = gsl_vector_get(x, i);
       _Complex double exp_iphi = cexp((-1.0*I)*phi);
      for(j=0;j<n;j++)
       {
           psi[j] += coef[i]*psiEs[j+base]*exp_iphi;
       }
    }
   return psi;
 }



//double L(_Complex double * arr, int n)
// {
//    double norm = 0.0;
//    int r;
//    for(r=0; r < n; r++)
//    norm += arr[r]*conj(arr[r]);
//  return norm;
// }

/*choosing states in which particles lie within x_init and x_fin*/
ull ** calc_regions_x(PARAMS * pm)

{
   int i;
   ull ** regions;
   ull * region_sizes;
   int n_reg = (pm->x_fin) - (pm->x_init) -1;
   ull * binary_basis = enumerate_r_basis(pm->num_sites,pm->num_particles);
   int num_in=0, num_out=0;
   int IN=0, OUT=1;

   newarr_(regions,2);

for(i=0; i < pm->numstates; i++)
   {
      unsigned long s = binary_basis[i];
      if (num_ones_in_range(pm->x_init, pm->x_fin, s) == pm->num_particles)
      {
         appendarr_(regions[IN],i);
      }
      else
      {
         appendarr_(regions[OUT],i);
      }
   }
   return regions;
}

/*providing a function, minex_func.f for the maximization process; calculating the probabilty using only localized states in calc_regions_x function*/
double my_prob_func(const gsl_vector * x, void * maxparams)
{
   PARAM_MAX * param_max = (PARAM_MAX *) maxparams; 
    double * psiEs;
    PARAMS * pm;
    CG * cg;
    double * evalues;
    ull ** regions;
    double prob = 0.0;
    double norm = 0.0;
    _Complex double * coef = param_max->coef;
    pm = param_max->pm;
    cg = param_max->cg;
    psiEs = param_max->psiEs;
    evalues = param_max->evalues;
    regions =param_max->regions;
    int j;

    _Complex double * psi = psi_phi_x(x, coef, psiEs);
    ull * region = regions[0]; 
    for(j=0; j < size_(region); j++)
      {
      // prob +=psi[region[j]]*conj(psi[region[j]]);  
       double * state = (double*) psi;  
       prob += SQR(state[2*region[j]])+SQR(state[2*region[j]+1]);
      }
   if (prob == 0)
      return 0;
   // return -prob/L(psi, n);     
return -prob;
}

/*Maximization process start here*/
int regional_prob_max(PARAMS * pm, CG * cg, _Complex double * coef ,double * psiEs, double * evalues, ull ** reg)
{
    const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex2;
    gsl_multimin_fminimizer *s = NULL;
    gsl_vector *ss, *x;
    gsl_multimin_function minex_func;
    int iter = 0;
    int status;
    double size;


    PARAM_MAX param_max;

    param_max.coef = coef;
    param_max.psiEs = psiEs;
    param_max.evalues = evalues;
    param_max.pm = pm;
    param_max.cg = cg;
    param_max.regions = reg;
    int n = size_(coef);
    int j;
    x = gsl_vector_alloc (n);

/*initializing the phases equally*/ 
   //gsl_vector_set_all (x, 0.0);

/*initializing the phases randomly*/
   for (j=0; j<n; j++)
    {
      double rand = ran() * 2 * M_PI;
      gsl_vector_set (x, j,rand);
    }

/* Set initial step sizes to 1 */
    ss = gsl_vector_alloc (n);
    gsl_vector_set_all (ss, 1.0);

/* Initialize method, iterate, and function to minimize */
    minex_func.n = n;
    minex_func.f = my_prob_func;
    minex_func.params = &param_max;
    s = gsl_multimin_fminimizer_alloc (T, n);
    gsl_multimin_fminimizer_set (s, &minex_func, x, ss);
    do
       {
         iter++;
         status = gsl_multimin_fminimizer_iterate(s);

         if (status)
             break;
         size = gsl_multimin_fminimizer_size (s);
         status = gsl_multimin_test_size (size, 1e-2);
         if (status == GSL_SUCCESS)
           {
             printf ("P_region converged to max at\n");
           }
         //printf ("%5d %10.3e f() = %7.3f size = %.3f\n",iter, gsl_vector_get(s->x, 0), s->fval, size);
       }

    while (status == GSL_CONTINUE && iter <2800000);
    printf ("%5d %10.3e f() = %7.3f size = %.3f\n",iter, gsl_vector_get(s->x,3) ,-1*(s->fval), size);

//calculate psi using the phase from the last iteration in the minimization loop (s->x)//
    _Complex double * psi = psi_phi_x(s->x, coef, psiEs);//calculate psi using the x from the last iteration in the minimization loop 
/*
//calculating entropy_EX given psi(s->x)//
    double S_o = ObsEntropyEX(pm, cg, psiEs, evalues, psi);
    printf ("S_EX from max_prob = %5f\n",S_o);
*/
//calculating entropy_ent given psi(s->x)//
    double Ent_ent = calc_ent_entropy_one_ev_complex_(psi, pm, pm->num_bath_sites);
    printf ("S_ent from max_prob = %5f\n",Ent_ent);
  
//calculating number density operator using psi(s->x) from the last iteration in the minimization loop//
    ull * binary_basis = enumerate_r_basis(pm->num_sites,pm->num_particles);
    double * dens_S_EX = den(pm, psi, binary_basis);
    printf("density_S_EX:\n");
    int index;
     for (index=0;index < pm->L;index++)
     {
        printf("%lf\n",dens_S_EX[index]); 
     }

     gsl_vector_free(x);
     gsl_vector_free(ss);
     gsl_multimin_fminimizer_free (s);
 
   return status;
 }

