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
      return -S_o;
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
     //printf("%lf\n",S_entang); 
     return -S_entang;
  
  }


//Minimization process start here//
double Entropy_min(PARAMS * pm, CG * cg, _Complex double *coef, double * psiEs, double * evalues,EG * eg, _Complex double * psi_init)
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
      double rr = ran();
      printf("randaom is %lf:",rr);
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
         status = gsl_multimin_test_size (size, 1e-2);
         
         if (status == GSL_SUCCESS)
           {
            printf ("S converged to minimum at\n");
           }   
         //printf ("%5d f() = %7.3f \n",iter, s->fval);
       }

    while (status == GSL_CONTINUE && iter < 8000000);
   double Smax=-(s->fval);
   printf ("Smax=%7.6f\n",Smax);
   //printf ("%7.7f\n",s->fval);
   //
   
 /*   
//outputting eigenvalues of reduced density matrix:
_Complex double * psi = psi_phi(s->x, coef, psiEs);
int L = pm->num_sites;
int np = pm->num_particles;
int bath = pm->num_bath_sites; 
unsigned  long NN = pm->numstates;
ull * bin = enumerate_r_basis(L,np);
int num_bath_particles;
int num_bath_particles_max;
num_bath_particles_max  = MIN(bath,np);
ull * basis_size =  init_bases(bath, L,np);


_Complex double ** reduced_rho = red_dens_mat_complex(L,np,bath, bin, NN, psi);

for(num_bath_particles=0; num_bath_particles < num_bath_particles_max+1; num_bath_particles++)
       {
         double * eigenvalues;
         double S_=0;
         int i, n = basis_size[num_bath_particles];
         calloc_(eigenvalues, n);
         diag_rho_cmplx(n, (double *) reduced_rho[num_bath_particles], eigenvalues);

         for(i=0; i < n; i++)
         {
            assert(eigenvalues[i] > -1.0e-10);
            if (eigenvalues[i] > 0)
            printf ("evalue[%f][%f]=%7.6f\n",eigenvalues[i],num_bath_particles,i);
         }
         free(eigenvalues);
        }

//Calculating the overlap btw |n,0> state and psi:
double p=0;
double phi_f;
int ii,jj;
int N = pm->numstates;
unsigned long long * binary_basis = enumerate_r_basis(pm->num_sites,pm->num_particles);

 for(jj=0;jj<N;jj++){ //going over all binary basis
  unsigned long b = binary_basis[jj];

  if (num_ones_in_range(0,pm->num_sites-pm->num_bath_sites , b) == pm->num_particles){
 //if (num_ones_in_range(pm->num_sites-pm->num_bath_sites,pm->num_sites, b) == pm->num_particles){
 
  _Complex double v = 0;
    for(ii=0; ii < N ;ii++){ //going over Evector
     double * evector = psiEs+N*ii;
     phi_f = gsl_vector_get(s->x, ii);//getting the phase for evector[ii] 
     _Complex double exp_iphi = cexp((1.0*I)*phi_f);
     v+=conj(coef[ii])*exp_iphi*evector[jj];}

  p += v*conj(v);}}
  printf("overlap with |n,0> is:%lf\n",p);

//printing all binary-states and probabilities:
double phi_f;
int ii,jj;
int N = pm->numstates;
unsigned long long * binary_basis = enumerate_r_basis(pm->num_sites,pm->num_particles);
double np_n0 = 0;
double np_0n = 0;
double np_nhalf = 0;
double np_12 = 0;
double np_21 = 0;
int sn0=0;
int s0n=0;
int sn2=0;
int s12=0;
int s21=0;
 for(jj=0;jj<N;jj++){ //going over all binary basis

  double p=0;
  unsigned long b = binary_basis[jj];
  _Complex double v = 0;
  for(ii=0; ii < N ;ii++){ //going over Evector
     double * evector = psiEs+N*ii;
     phi_f = gsl_vector_get(s->x, ii);//getting the phase for evector[ii] 
     _Complex double exp_iphi = cexp((1.0*I)*phi_f);
     v+=conj(coef[ii])*exp_iphi*evector[jj];}
    p += v*conj(v);
    //print_binary(b, pm->num_sites);
    //printf(" binary: %ld\n",b);
    //printf("P for binary %lu is = %lf\n",b,p);
    if (num_ones_in_range(0,pm->num_sites-pm->num_bath_sites,b) == pm->num_particles){
    np_n0+=p;
    sn0+=1;}
    if (num_ones_in_range(pm->num_sites-pm->num_bath_sites,pm->num_sites,b) == pm->num_particles){
    np_0n+=p;
    s0n+=1;}
    //if (num_ones_in_range(pm->num_sites-pm->num_bath_sites,pm->num_sites, b) == (pm->num_particles)/2){
    //np_nhalf+=p;
    //sn2+=1;}
    
    if (num_ones_in_range(pm->num_sites-pm->num_bath_sites,pm->num_sites,b) == 2){
    np_12+=p;
    s12+=1;}
    if (num_ones_in_range(pm->num_sites-pm->num_bath_sites,pm->num_sites,b) == 1){
    np_21+=p;
    s21+=1;}
    
    }
    printf("P|n,0>=%lf , n=%d\n",np_n0, sn0);
    printf("P|0,n>=%lf , n=%d\n",np_0n, s0n);
    //printf("P|n/2,n/2>=%lf,n=%d\n",np_nhalf, sn2);
    printf("P|1,2>=%lf , n=%d\n",np_12, s12);
    printf("P|2,1>=%lf , n=%d\n",np_21, s21);



//Calculating the coarse_density//
  _Complex double * psi = psi_phi(s->x, coef, psiEs); 
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
if ((s->fval)<-2.3){
 double np=0;
 //printf("S max is:%lf\n",Smax);
printf("S max is:%lf\n",s->fval); 
_Complex double * psi = psi_phi(s->x, coef, psiEs);
 ull * binary_basis = enumerate_r_basis(pm->num_sites,pm->num_particles);
 double * density_matrix = den(pm, psi, binary_basis);
 printf("density_S:\n");
 int index;
 for (index=0;index < pm->L;index++)
    {
       if (index >=0 && index<((pm->L)-(pm->num_bath_sites))){
         np += density_matrix[index];}
         //printf("%lf\n",density_matrix[index]);
    }
         printf("number of particles in the subsystem is: %lf\n", np);
         double np_bath=pm->num_particles - np;
         //printf("number of particles in bath is: %lf\n", np_bath);
         double D = np_bath-np;
         //printf("D = np_bath-np: %lf\n",D);
   }       
        
//calculating S_ent using s->x // 
//_Complex double * vec = psi_phi(s->x, coef, psiEs);
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


//calculating psif using psiEs and region

         int l,jj,f;
         complx * psif;
         newarr_(psif,N);
         ull * regions = (calc_regions(pm))[0];//chnage this to some other region where it can be controled by x_init and x_fin 
         int size = size_(regions);
         double phi_f;
         for (f=0;f<size;f++){
         for (l=0;l<N;l++){
         phi_f = gsl_vector_get(s->x, l);
         _Complex double exp_iphif = cexp((-1.0*I)*phi_f);
         double * evector = psiEs+N*l;
         psif[regions[f]] += coef[l]*exp_iphif*evector[regions[f]];}


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
/*
_Complex double ** reduced_dens_mat_complex(int num_sites, int num_particles, int 
num_bath_sites, ull * binary_basis, unsigned  long nbasis, _Complex double * evector)
{

   _Complex double ** rho;
   int num_bath_particles_max  = MIN(num_bath_sites, num_particles);
   int i;
   unsigned long mask;
   unsigned long base=0;
   ull * basis_size =  init_bases(num_bath_sites, num_sites, num_particles);

   rho = callocate_rho_cmplx(num_bath_particles_max, basis_size);
   mask  = ((1 << num_bath_sites) -1) << (num_sites-num_bath_sites);


   while (base < nbasis)
 {
      int nbath = num_ones(mask&binary_basis[base]);
      int count = 0;
      if (basis_size[nbath] == 0)
        {
          base++;
          continue;
        }
      for(i=0; i < basis_size[nbath];i++)
      {
         int j;
         for(j=0; j <= i;j++)
         {
            rho[nbath][count] += evector[base+i]*conj(evector[base+j]);
            count++;
         }
      }
      base += basis_size[nbath];
   }
   free(basis_size);
   return rho;}*/
