#include "locate_multivs.h"
//#include <gsl/gsl_multimin.h>
//#include <gsl/gsl_vector.h>
#include "min.h"
#include "ent_entropy.h"
#include "energy_graining.h"
#include "monte_carlo.h"
#include <time.h>
#include <stdio.h>
#include "defs.h"
int * distinct_rands(int n, int range_l, int range_h)
{
   int * mark=0, * out_list=0;
   int num=0;
   newarr_(mark,range_h);
   while (num < n)
   {
      int new_rand = randint(range_l,range_h);
      if (!mark[new_rand])
      {
         mark[new_rand] = 1;
         appendarr_(out_list,new_rand);
         num++;
      }
   }
   freearr_(mark);
   qsort(out_list,size_(out_list),sizeof(out_list[0]),int_comp);
   return out_list;
}

_Complex double * psi_phi_change(int index, double phi_new, double * phi, _Complex double * coef, double * psiEs, _Complex double * psi, int change)
{
   int n = size_(coef);
   int base = n*index;
   _Complex double Delta_exp_iphi = cexp((-1.0*I)*phi_new)-cexp((-1.0*I)*phi[index])*change;

   int j;
   for(j=0;j<n;j++)
   {
      psi[j] += coef[index]*psiEs[j+base]*Delta_exp_iphi;
   }
   return psi;
}

double S_ent(_Complex double * psi, void * params)
{
   PARAM_MIN * param_min = (PARAM_MIN *) params;

   PARAMS * pm = param_min->pm;
   CG * cg = param_min->cg;
   _Complex double * coef = param_min->coef;
   double * psiEs = param_min->psiEs;
   double S_entang = calc_ent_entropy_one_ev_complex_(psi, pm, pm->num_bath_sites);
   return S_entang;

}
double S_xe(_Complex double * psi, void * params)
{  
   PARAM_MIN * param_min = (PARAM_MIN *) params;
   
   PARAMS * pm = param_min->pm;
   CG * cg = param_min->cg; 
   _Complex double * coef = param_min->coef;
   double * psiEs = param_min->psiEs;
   double * evalues = param_min->evalues;
   double S_obs = ObsEntropyEX(pm, cg, psiEs, evalues, psi);
   free2darr_(cg->Ps);
   free2darr_(cg->density);
   return S_obs;

}

double monte_func(PARAMS * pm, CG * cg, _Complex double *coef, double * psiEs, double * evalues,EG * eg, _Complex double * psi_init)
{
    clock_t start = clock();
    FILE * monte_carlo_entropy_data = fopen("monte_carlo_entropy_data.d","w");
    double  * phi;
    int iteration;
    double c =0;
    int accept = 0;
    PARAM_MIN param_min;
    param_min.coef = coef;
    param_min.psiEs = psiEs;
    param_min.evalues = evalues;
    param_min.pm = pm;
    param_min.cg = cg;
    param_min.eg = eg;
    int n = size_(coef);
    int bath = pm->num_bath_sites;
    int iter = pm->iter;
    double step_size = pm->step_size;
    double beta_carlo = pm->beta_carlo;
    _Complex double  * psi;
    
    newarr_(phi, n);
    newarr_(psi, n);

    
    //Initializing the phases randomly//
    int i,j;
    for (j=0; j<n; j++)
    {
        double rand = ran() * 2 * M_PI;
        phi[j] = rand;
    }
    
    for (j=0; j<n; j++)
    {
        psi_phi_change(j, phi[j],  phi, coef, psiEs, psi, 0);
    }
    printf("in monte_func: initially  L2(psi) = %lf\n",L2((double *) psi,2*n));
   
    double S_init;
    if (pm->monte_ent == 1)   S_init = S_ent(psi, &param_min);
    if (pm->monte_Sxe == 1)   S_init = S_xe(psi, &param_min);
    if (pm->monte_density == 1){
    ull * binary_basis = enumerate_r_basis(pm->num_sites,pm->num_particles);
    double * density_init = den(pm, psi, binary_basis);
    S_init = density_init[0];}

    _Complex double * psi_old;
    newarr_(psi_old,n);
    int * indices=0;
    double * phi_new, * phi_inits;
    int m=pm->monte_moves;//= number of phases changed in one step
    newarr_(phi_new,m);
    newarr_(phi_inits,m);

    for (iteration=0; iteration < iter; iteration++){
       int j;
       //pick random index
       indices = distinct_rands(m,0,n);

       memcpy_(psi_old,psi,n);
       //saving the old phi[i] in phi_init
       for(j=0;j<m;j++) {
	  i = indices[j];
	  phi_inits[j]= phi[i];
	  double delta_phi =  2*M_PI*(2*ran()-1)*step_size;
	  //randomly change phi_i
	  phi_new[j] =  fmod(phi_inits[j]+delta_phi,(2*M_PI));
	  psi = psi_phi_change(i, phi_new[j],  phi, coef, psiEs, psi, 1);
	  phi[i] =  phi_new[j];
       }

       double S_attempt; 
       if (pm->monte_ent == 1)  S_attempt = S_ent(psi, &param_min);
       if (pm->monte_Sxe == 1)  S_attempt = S_xe(psi, &param_min);
       if (pm->monte_density == 1){
	  ull * binary_basis = enumerate_r_basis(pm->num_sites,pm->num_particles);
	  double * density_init = den(pm, psi, binary_basis);
	  S_attempt = density_init[0];}
       double deltaS = S_attempt - S_init;

       //condition for rejecting
       double dd = ran();
       if (exp(- beta_carlo * deltaS) < dd){ //reject
	  for(j=0;j<m;j++) {
	     i = indices[j];
	     phi[i] = phi_inits[j];		
	  }
	  memcpy_(psi,psi_old,n);
       }
       else {
	  accept += 1;
	  S_init = S_attempt;
       }
       
       if (1 || iteration % 1000 ==0) {
         fprintf(monte_carlo_entropy_data,"%lf\n", S_init);
	//fprintf(monte_carlo_entropy_data,"%d %lf\n",iteration ,S_init);

       }

       if (iteration % 10000 ==0) {
	  c = (float)accept/(float)(iteration+1);
	  printf("acceptance ratio at iteration %d is %lf\n",iteration, c);
          printf("in monte_func: L2(psi) = %lf\n",L2((double *) psi,2*n));
	  //printf("%d %lf\n",iteration, c);
       }
       freearr_(indices);
    }
    freearr_(phi_new);
    freearr_(phi_inits);
    c = (float)accept/(float)(iteration+1);
    fclose(monte_carlo_entropy_data);
    freearr_(psi);
    freearr_(phi);
    freearr_(psi_old);
    clock_t stop = clock();
    double elapsed = (double)(stop - start) * 1000.0 / CLOCKS_PER_SEC;
    printf("Time elapsed in ms: %f\n", elapsed);
    return c;  
}	



















