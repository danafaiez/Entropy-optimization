//This program minimized factorized observational entropy, observational entropy with position and energy coarse graining, or entanglement entropy//
#include "locate_multivs.h"
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_vector.h>
#include "min.h"
#include "ent_entropy.h"
#include "energy_graining.h"
//Unitary evolution of the wavefunction//
_Complex double * psi_phi(const gsl_vector *x, _Complex double * coef, double * psiEs)
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
//providing a function, minex_func.f for the minimization process; calculating factorized observational entropy//
double my_f_fac(const gsl_vector * x, void * params)
{
     PARAM_MIN * param_min = (PARAM_MIN *) params;
     complex double * coef;
     PARAMS * pm;
     CG * cg;
     double * psiEs;
     double * evalues;
     EG * eg;
     coef = param_min->coef;
     pm = param_min->pm;
     cg = param_min->cg;
     psiEs = param_min->psiEs;
     evalues = param_min->evalues;
     eg = param_min->eg;
    _Complex double * psi = psi_phi(x, coef, psiEs); 
    _Complex double * psi_e_b = transform_pos_to_energy(eg, psi);
    double S_f = Sobs_fine_grain_E(psi_e_b);
    freearr_(psi);
    return S_f;
}
//providing a function, minex_func.f for the minimization process; calculating observational entropy with position and energy coarse graining//
double my_f_ex(const gsl_vector * x, void * params)
{
    PARAM_MIN * param_min = (PARAM_MIN *) params;
 
    _Complex double * coef;
    PARAMS * pm;
    CG * cg; 
    double * psiEs;
    double * evalues;
    coef = param_min->coef;
    pm = param_min->pm; 
    cg = param_min->cg; 
    psiEs = param_min->psiEs;
    evalues = param_min->evalues;

      _Complex double * psi = psi_phi(x, coef, psiEs);
      double S_o = ObsEntropyEX(pm, cg, psiEs, evalues, psi);
      free2darr_(cg->Ps);
      free2darr_(cg->density);
      freearr_(psi);
      return S_o;
}
//providing a function, minex_func.f for the minimization process; calculating entanglement entropy//
double my_f_ent(const gsl_vector * x, void * params)
  {
      PARAM_MIN * param_min = (PARAM_MIN *) params;
  
      _Complex double * coef;
      PARAMS * pm;
      CG * cg;
      double * psiEs;
      double * evalues;
      coef = param_min->coef;
      pm = param_min->pm;
      cg = param_min->cg;
      psiEs = param_min->psiEs;
      evalues = param_min->evalues;
      _Complex double * psi = psi_phi(x, coef, psiEs);
      double S_entang = calc_ent_entropy_one_ev_complex_(psi, pm, pm->num_bath_sites);
      freearr_(psi);
  
     return S_entang;
  
  }


//Minimization process start here//
double Entropy_min(PARAMS * pm, CG * cg, _Complex double *coef  ,double * psiEs, double * evalues,EG * eg)
{
    const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex2;
    gsl_multimin_fminimizer *s = NULL;
    gsl_vector *ss, *x;
    gsl_multimin_function minex_func;
    int iter = 0;
    int status;
    double size;


    PARAM_MIN param_min;   
    param_min.coef = coef;
    param_min.psiEs = psiEs;
    param_min.evalues = evalues;
    param_min.pm = pm;
    param_min.cg = cg;
    param_min.eg = eg;
    int n = size_(coef);
  
     x = gsl_vector_alloc (n);
//Initializing the phases equally//
//   gsl_vector_set_all (x, 0.0);

//Initializing the phases randomly//
     int j;
     for (j=0; j<n; j++)
     {
      double rand = ran() * 2 * M_PI;
      gsl_vector_set (x, j,rand);
     }
//Set initial step sizes to 1 //
    ss = gsl_vector_alloc (n);
    gsl_vector_set_all (ss, 0.2);

//Initialize method, iterate, and function to minimize//
    minex_func.n = n;
    if (pm->minimize_S_EX == 1)  minex_func.f = my_f_ex;
    
    if (pm->minimize_Sent == 1)   minex_func.f = my_f_ent;
 
    if (pm->minimize_foe == 1)   minex_func.f = my_f_fac;
    
    minex_func.params = &param_min;  
    s = gsl_multimin_fminimizer_alloc (T, n);
    gsl_multimin_fminimizer_set (s, &minex_func, x, ss);
    do
       {
         iter++;
         status = gsl_multimin_fminimizer_iterate(s);
        
         if (status)
             break;
         size = gsl_multimin_fminimizer_size (s);
         status = gsl_multimin_test_size (size, 1e-9);
         
         if (status == GSL_SUCCESS)
           {
            //printf ("S converged to minimum at\n");
           }   
         //printf ("%5d f() = %7.3f \n",iter, s->fval);
       }

    while (status == GSL_CONTINUE && iter < 2000000);
    // printf ("Smin = %7.6f\n",s->fval);
      printf ("%7.7f\n",s->fval);//}


//Calculating the coarse_density//
//calculate psi using the phase from the last iteration in the minimization loop (s->x)//
/*  _Complex double * psi = psi_phi(s->x, coef, psiEs); 
    double ** Ps = calc_PsEX(pm, cg, evalues, psiEs, psi);
    double ** density = coarse_density(pm, cg, Ps);
    double sum=0;
    int n_boxes=cg->L/cg->size_of_box;
    int E_index,j;
    int ncoarse = n;
    
    for (E_index=0;E_index < ncoarse;E_index++)
   {
      for(j=0; j < n_boxes;j++)
      {
       printf(" %lf " , density[E_index][j]); 
       sum += density[E_index][j];
      }
    printf("\n");  
    
   }
printing the sum of the probabilities*
   printf("%lf", sum);
   printf("\n");


//calculating number density operator using psi(s->x) from the last iteration in the minimization loop//

//condition for calculating <N>//
//if ((s->fval)<0.007){

 _Complex double * psi = psi_phi(s->x, coef, psiEs);
 ull * binary_basis = enumerate_r_basis(pm->num_sites,pm->num_particles);
 double * density_matrix = den(pm, psi, binary_basis);
 printf("density_S:\n");
 int index;
 for (index=0;index < pm->L;index++)
    {
        printf("%lf\n",density_matrix[index]);
    }
//}

//calculating S_ent using s->x // 
_Complex double * vec = psi_phi(s->x, coef, psiEs);
//double S_ent_corres = calc_ent_entropy_one_ev_complex_(vec, pm, pm->num_bath_sites);

//calculating S_ex using s->x // 
_Complex double * vec = psi_phi(s->x, coef, psiEs);
double S_ex_corres = ObsEntropyEX(pm, cg, psiEs, evalues, vec);
printf("S_ex_corres = %lf\n",S_ex_corres);

//calculating FOEusing s->x // 
_Complex double * psi_e_b_corres = transform_pos_to_energy(eg, vec);
 double S_f_corres = Sobs_fine_grain_E(psi_e_b_corres);

printf("S_ent_corres = %lf\n",S_ent_corres);
printf("S_ex_corres = %lf\n",S_ex_corres);
printf("S_f_corres = %lf\n",S_f_corres);


//calculating expectation value of energy for the state with minS
//calculating psi1 using psiEs and region
//if((s->fval)<2.75){
         int i,u;
         double p,expE=0;
         complx * psi1;
         newarr_(psi1,n);
        ull * regions = (calc_regions(pm))[0];
        int sz = size_(regions);
        for (i=0;i<n;i++){
        u=0;
          while (u<sz)
          if (i==regions[u]){
          for (j=0;j<n;j++){
          p = gsl_vector_get(s->x, j);
          _Complex double exp_ip = cexp((-1.0*I)*p);
          double * evector = psiEs+n*j;
          psi1[i] += coef[j]*evector[i]*exp_ip;}
          break;}
          else{
           u++;
           psi1[i]=0;}}

          int ii;
          _Complex double * c;
          double exp_E=0;
          newarr_(c,n);
          c = coeff(psi1, size_(psi1), psiEs);
          for(ii=0;ii<n;ii++){
          exp_E += evalues[ii]*SQR(cabs(c[ii]));}
          printf("Energy of small box:\n");
          printf("%lf\n",exp_E);
//}

*/
    gsl_vector_free(x);
    gsl_vector_free(ss);
    gsl_multimin_fminimizer_free(s);
//  return status;
    double sfval = s->fval;
    return sfval;
}
//calculating <N>// 

  double * den(PARAMS * pm, _Complex double * psi, ull *  binary_basis)
   {
       int i;
       double * rho;
       newarr_(rho,pm->L);
       for(i=0;i < pm->numstates;i++)
       {
            int j;
            double re=creal(psi[i]);
            double im=cimag(psi[i]);
            double prob = re*re+im*im;
 
       for(j=0; j < pm->L; j++)
         {
               ull occupied =  bit_val(binary_basis[i],j);
               rho[j] += occupied*prob;
         }
      }
  return rho;
  }



