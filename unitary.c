#include "locate_multivs.h"
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_multimin.h>
#include "energy_graining.h"
#include "basis_states.h"
#include "unitary.h"
#include "ran.h"
#include "min.h"
#define THETA 0
#define PHI 1
#define PSI 2

/*
typedef struct {
   double theta,phi,psi;
} ANGLES;
*/

#define Mflat(i,j)  Mat[i*N+j]

double * packMat (double * Mat,int N)
{
   double * Mpacked;
   newarr_(Mpacked,(N*(N+1))/2);
   int i,j,indx=0;
   for(i=0;i < N;i++)
      for(j=0; j <= i;j++)
         Mpacked[indx++] = Mflat(j,i);
   return Mpacked;
}

double * firstDeriv(complx *s, double **J)
{
   int N = size_(J);
   int i,j;
   double * derivs;
   newarr_(derivs,N);
   for(i=0;i<N;i++)
   {
      complx d=0;
      for(j=0;j<N;j++)
      {
         d += J[i][j]*s[j];
      }
      derivs[i] = -2*cimag(conj(s[i])*d);
   }
   return derivs;
}

int diag_rho(int N, double * matrix_double, double * eigenvalues);

double *  scndDerivEV(complx * s,double ** J)
{

   double * Mat;
   double * Mpacked;
   int N = size_(J);
   newarr_(Mat,N*N);
   int i,j;
   printf("***making 2nd deriv matrix...\n");
   for(i=0;i < N; i++)
   {
      double dia=0;
      for(j=0;j < N; j++)
      {
         Mflat(i,j) = creal(s[i]*conj(s[j]))*J[i][j];
         dia += Mflat(i,j);
      }
      Mflat(i,i) -= dia;
   }

   for(i=0;i < N; i++)
     for(j=0;j < N; j++)
        Mflat(i,j) = -Mflat(i,j);

   Mpacked = packMat(Mat,N);
   

   printf( "diagonalizing...\n");
   double * eigenvalues;
   newarr_(eigenvalues,N);
   diag_rho(N, Mpacked, eigenvalues);
   freearr_(Mat);
   freearr_(Mpacked);
   return eigenvalues;
}

#undef Mflat


complx * compute_partial_fields(double * a, complx * s,int N,int M,complx ** E)
{
   complx * h;
   newarr_(h,M);
   int j,k;
   for(k=0;k < M;k++)
      for(j=0;j < N;j++)
      {
         h[k] += a[j]*E[j][k]*conj(s[j]);
      }

   return h;
}


int test_2nd_deriv(PSI_STATE * psi_state, complx * s2)
{
   int minimum_found = 0;
   int N = psi_state->N;
   int M = psi_state->M;
   int k;
   complx ** E = psi_state->E;
   complx ** newE = psi_state->newE;
   double * a = psi_state->a;
   double * evs=0;

   complx * h =  compute_partial_fields(a, s2,N,M,newE);
   double err=0;
   for(k=1;k < M;k++)
      err += h[k]*conj(h[k]);
   err = sqrt(err);
   printf("error in partial fields = %lf\n", err);
   if (err > 0.1)
      minimum_found = -1;
   if (err < 1e-5)
   {

      evs = scndDerivEV(s2,psi_state->J);
      if (evs[0] > -1e-10)
      {
         minimum_found = 1;
      }
      freearr_(evs);
   }
   freearr_(h);
   return minimum_found;
}


double  newPsi(PSI_STATE * psi_state)
{
   
   int numstates = psi_state->numstates; 
   int size_of_box = psi_state->size_of_box;
   PARAMS * pm = psi_state->pm;
   double * psiEs = psi_state->psiEs;
   double * energy = psi_state->energy;
   CG * cg = psi_state->cg;
  
   int i,j,count,countmax=40000;
   int N = psi_state->N;
   int M = psi_state->M;
   complx *z;
   complx ** E = psi_state->E;
   complx ** newE = psi_state->newE;
   double * a = psi_state->a;
   complx * psi;
   complx * psiN;
   double norm;
   newarr_(psi,M); 
   newarr_(z,N);

   for(i=0;i<N;i++)
   {
      complx g = c_gaussian_random();
      z[i] = g/cabs(g);
   }

   for(count=0;count < countmax;count++)
   {
      norm = 0;
      int x,iE;
      for(x=0;x < M;x++)
      {
         psi[x]=0;
         for(iE=0;iE<N;iE++)
         {
            psi[x] += a[iE]*z[iE]*E[iE][x];
         } 
         norm += SQR(creal(psi[x]))+SQR(cimag(psi[x]));
      }
      if (count%200 == 0)
      {
         int minimum_found=0;
         printf("P = %lf\n",norm);
         double * derivs = firstDeriv(z,psi_state->J);
         double deriv_error = L2(derivs,N);
         deriv_error = sqrt(deriv_error);
         printf("deriv_error = %g\n",deriv_error);
         
         if (deriv_error < 1e-10)
         {
           minimum_found=1;
           if (0)
             {   
            double * evs = scndDerivEV(z,psi_state->J);
            if (evs[0] > -1e-10)
            {
               printf("minimum found\n");
               minimum_found = 1;
            }
            freearr_(evs);
            }
         }
        
        if (minimum_found){
       //observational entropy//       
         newarr_(psiN,N);
         int y;      
        for(y=0;y < N;y++)
        {
         for(iE=0;iE<N;iE++)
         {
          double * evector = psiEs+N*iE; 
          psiN[y] += a[iE]*z[iE]*evector[y];
         }
       }
         //CG * cg = create_CG(pm, size_of_box,numstates);
         double S_o = ObsEntropyEX(pm, cg, psiEs, energy, psiN);
         printf ("S_EX(unitary maxP) = %5f\n",S_o); 

//calculating number density
      
         ull * binary_basis = enumerate_r_basis(pm->num_sites,pm->num_particles);
         double * density_matrix = den(pm, psiN, binary_basis); 
         printf("density_S:\n");
         int index;
         for (index=0;index < pm->L;index++){
         printf("%lf\n",density_matrix[index]);}
     
//calculating psi1 using psiEs and region
/*
         int i,j,u,J,l,ii;
         double expE=0;
	 complx * psi1;
         newarr_(psi1,N);
 
	ull * regions = (calc_regions(pm))[0];
	int size = size_(regions);
	for (i=0;i<N;i++){
	u=0;
        while (u<size)
	if (i ==regions[u]){
	for (j=0;j<N;j++){
	double * evector = psiEs+N*j;
	psi1[i] += a[j]*z[j]*evector[i];}  break;}  //j loop closed; if loop closed 
	else{
	if (u==size){
	psi1[i]=0;
	break;}
	else
	u++;}} //else loop closed, i loop closed


*/
//calculating psi1 using W  
  
        double expE=0; 
        int u,J,l,ii;
        complx * psi1;
        newarr_(psi1,N);
        complx ** W = makeEN(pm, psiEs);
  	for(u=0;u<N;u++){
        for(J=0;J<N;J++){
       psi1[u] += a[J]*z[J]*W[J][u];}}        
          

//calculating <E>  
          _Complex double * c;
          newarr_(c,N);
          c = coeff(psi1, size_(psi1), psiEs);
          for(ii=0;ii<N;ii++){
          expE += energy[ii]*SQR(creal(c[ii]))+SQR(cimag(c[ii]));}
          printf("Energy of small box:\n");
          printf("%lf\n",expE);
          freearr_(derivs);

//calculating entanglement entropy....


          return norm;}
      }
      norm = 1.0/sqrt(norm);
      for(x=0;x < M;x++)
         psi[x] *= norm;

      for(iE=0;iE<N;iE++)
      {
	 z[iE]=0;
         for(x=0;x<M;x++)
         {
            z[iE] += conj(E[iE][x])*psi[x];
         }
         z[iE] = z[iE]/cabs(z[iE]);
      }
   }
   //freearr_rr_(psi);
   return 1e10;
}

complx ** makeEN(PARAMS * pm, double * evectors)
{
   unsigned long long * binary_basis = enumerate_r_basis(pm->num_sites, pm->num_particles);
   ull * regions = (calc_regions(pm))[0];
   int i,j;
   int N = pm->numstates;
   complx ** W;
   new2darr_(W,N,N);
   int x_begin= (N - pm->num_bath_sites)/2;
   for(i=0;i<N;i++)
   {
   double * evector = evectors+N*i;
   for(j=0;j<N;j++)
   {
      unsigned long s = binary_basis[j];
      if (num_ones_in_range(0, pm->num_bath_sites, s) == pm->num_particles)
   //   if (num_ones_in_range(x_begin, x_begin+pm->num_bath_sites, s) == pm->num_particles)
       {
      W[i][j] = evector[j];
   //   W[i][j] = 0;
        }
      else
       {
     W[i][j] = 0;   
     // W[i][j] = evector[j];
       } 
   }
   }
   return W;
}

/*
complx ** makeEN(PARAMS * pm, double * evectors)
{
  // unsigned long long * binary_basis = enumerate_r_basis(pm->num_sites, pm->num_particles
   ull * regions = (calc_regions(pm))[0];
   int i,j,u;
   int N = pm->numstates;
   complx ** W;
   new2darr_(W,N,N);
   int x_begin= (N - pm->num_bath_sites)/2;
   for(i=0;i<N;i++)
   {
   double * evector = evectors+N*i;
   for(j=0;j<N;j++)
   {
      //for(u=u;u<size_(regions);u++){ 
       if (regions[]==j)
        {
      W[i][j] = evector[j];
   //   W[i][j] = 0;
        }
      else
       {
     W[i][j] = 0;   
     // W[i][j] = evector[j];
       } 
   }
   }
   return W;
}

*/
double ** makeJ(PSI_STATE * psi_state)
{
   int i,j,k;
   int N = psi_state->N;
   int M = psi_state->M;
   double ** J;
   complx ** E = psi_state->E;
   double * a = psi_state->a;
   new2darr_(J,N,N);

   for(i=0;i < N; i++)
     for(j=0;j < N; j++)
       for(k=0;k < M; k++)
          J[i][j] += creal(a[i]*a[j]*E[i][k]*conj(E[j][k]));
   return J;
}

complx * make_unitary_matrix2(double * angles)
      //double theta,double phi,psi
{
      complx * arr;
      complx exp_phi = cexp(angles[PHI]*1j);
      complx exp_phiC = conj(exp_phi);
      complx exp_psi = cexp(angles[PSI]*1j);
      complx exp_psiC = conj(exp_psi);
      double c = cos(angles[THETA]);
      double s = sin(angles[THETA]);
      newarr_(arr,4);
      arr[0] =  c*exp_phi;
      arr[1] = s*exp_psi;
      arr[2] = -s*exp_psiC;
      arr[3] = c*exp_phiC;
      return arr;
}

complx ** make_trans_state2(PSI_STATE * psi_state, double * angles ,int i, int
j)
{
   complx ** newE = psi_state->newE;
   complx ** E = psi_state->E;
   int k;
   int N = psi_state->N;
   int M = psi_state->M;
   complx * arr;
   arr = make_unitary_matrix2(angles);
   put2darr_(newE,E); 
   for(k=0; k < N; k++)
   {
      newE[k][i] = arr[0]*E[k][i]+ arr[1]*E[k][j];
      newE[k][j] = arr[2]*E[k][i]+ arr[3]*E[k][j];
   }
   freearr_(arr);
   return newE;
}

double energy_complex(PSI_STATE * psi_state, complx * s)
{
   int N = psi_state->N;
   int M = psi_state->M;
   double * a = psi_state->a;
   complx ** newE = psi_state->newE;
   complx ** E = psi_state->E;
   int i,k;
   double e=0;
   for(k=0; k < M; k++)
   {
      complx e_site=0;
      for(i=0; i < N; i++)
      {
         e_site += a[i]*conj(s[i])*newE[i][k];
         //e_site += a[i]*conj(s[i])*E[i][k];
      }
      e += SQR(creal(e_site))+SQR(cimag(e_site));
   }
   return -e;
}
double energy_complex1(PSI_STATE * psi_state, complx * s)
{
   int N = psi_state->N;
   int i,j;
   complx e=0;
   double ** J = psi_state->J;
   for(i=0; i < N; i++)
      for(j=0; j < N; j++)
         e += s[i]*J[i][j]*conj(s[j]);
   return -creal(e);
}



double energy_complex2(PSI_STATE * psi_state, complx * s)
{
   int N = psi_state->N;
   int M = psi_state->M;
   double * a = psi_state->a;
   complx ** newE = psi_state->newE;
   complx ** E = psi_state->E;
   int i,k;
   double e=0;
   double e_new=0;
   for(k=0; k < M; k++)
   {
      complx e_site=0;
      complx e_site_new=0;
      for(i=0; i < N; i++)
      {
         /*
         e_site += a[i]*conj(s[i])*E[i][k];
         e_site_new += a[i]*conj(s[i])*newE[i][k];
         */

         e_site += a[i]*conj(s[i]*s[i])*E[i][k];
         e_site_new += a[i]*conj(s[i]*s[i])*newE[i][k];
      }
      e += SQR(creal(e_site))+SQR(cimag(e_site));
      e_new += SQR(creal(e_site_new))+SQR(cimag(e_site_new));
   }
   if (fabs(e-e_new) > 1e-8)
   {
      printf("debugging energy_complex: e = %lf, e_new = %lf\n",e,e_new);
   }
   return -e;
}

complx * normali(complx * s)
{
   complx * s0;
   int N = size_(s0);
   int k;
   newarr_(s0,N);
   for(k=0; k < N;k++)
      s0[k] = s[k]/cabs(s[k]);
   return s0;
}

complx * normali2(complx ** E, int i)
{
   complx * s0;
   int N = size_(E);
   int k;
   newarr_(s0,N);
   for(k=0; k < N;k++)
      s0[k] = E[k][i]/cabs(E[k][i]);
   return s0;
}

//double minEunitary(PSI_STATE * psi_state, ANGLES * angles,i,j)
double minEunitary( const gsl_vector * x, void * params)
{

   PARAM_UNITARY * param_unitary = (PARAM_UNITARY *) params;
   PSI_STATE * psi_state = param_unitary->psi_state;
   int i = param_unitary->i;
   int j = param_unitary->j;
   int N = psi_state->N;
   int M = psi_state->M;
   double ** J = psi_state->J;
   int k;
   double angles[3];
   for(k=0;k<3;k++)
      angles[k] = gsl_vector_get(x,k);

   complx ** newE = make_trans_state2(psi_state,angles,i,j);
   complx * s0, s0C;
   s0 = normali2(newE,0);
   s0C = conj(s0[0]);
   for(k=0; k < N;k++)
      s0[k] *= s0C;
   double e = energy_complex(psi_state,s0);
   double e1 = energy_complex1(psi_state,s0);
   if (0 && fabs(e-e1) > 1e-8)
   {
      printf("2 energy methods differ, %lf, %lf\n",e,e1);
   }
   //   free2darr_(newE);
   freearr_(s0);
   return e;
}

double make_ran_state(int N,int M, complx ** E,complx ** newE, double * a,
double **J)
{
   int num_descents = 0;
   int max_num_descents = 20;
   double * evs=0;
   PSI_STATE psi_state;

   PARAM_UNITARY param_unitary;
   psi_state.N = N;
   psi_state.M = M;
   psi_state.E = E;
   psi_state.newE = newE;
   psi_state.a = a;
   psi_state.J = J;

   param_unitary.psi_state = &psi_state;

   const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex2;
   gsl_multimin_fminimizer *s = NULL;
   gsl_vector * ss, *x;//ss=stepsizes
   gsl_multimin_function minex_func;
   int i,j;
   int iter = 0;
   int status;
   x = gsl_vector_alloc (3);
   minex_func.n = 3;
   minex_func.params = &param_unitary;  
   minex_func.f = minEunitary;
   minex_func.params = &param_unitary;  

   //Set initial step sizes to 1 //
   ss = gsl_vector_alloc (3);

   while (num_descents < max_num_descents)
   {
      gsl_vector_set_all (ss, 1.0);
      s = gsl_multimin_fminimizer_alloc (T, 3);
      num_descents += 1;
      int count = 0;
      int countmax = 40000;
      double min_so_far = 1e10;
      int ij=0;
      while (count < countmax)
      {
         double p0[3];
         count += 1;
         //i = randint(0,M);
         i = 0;
         while (1)
         {
            j = randint(0,M);
            if (i != j)
               break;
         }
         param_unitary.i = i;
         param_unitary.j = j;
         int k;
         for(k=0;k < 3;k++)
         {
            p0[k] = 1.5*ran();
            gsl_vector_set (x, k, p0[k]);
         }
         gsl_multimin_fminimizer_set (s, &minex_func, x, ss);
         iter=0;
         printf("count = %d, i=%d j=%d\n",count,i,j);
         do
         {
            iter++;
            //printf("iter = %d\n",iter);
            status = gsl_multimin_fminimizer_iterate(s);

            if (i*j !=0 || status)
               break;
            double size = gsl_multimin_fminimizer_size (s);
            status = gsl_multimin_test_size (size, 1e-9);

            if (status == GSL_SUCCESS)
            {
               printf ("S converged to minimum %g, min_so_far = %g\n",s->fval,
min_so_far);
            }   
            //printf ("%d [%g,%g,%g] f() = %g size =
            //%g\n",iter,gsl_vector_get(s->x,0),gsl_vector_get(s->x,1),gsl_vector_get
            //(s->x,2),s->fval,size);
         }
         while (status == GSL_CONTINUE && iter < 1000);


         if (s->fval < min_so_far)
         {
            min_so_far = s->fval;
            printf ("min_so_far %g\n",min_so_far);
            put2darr_(E,newE);
         }
         if (count%100 == 0)
            printf("count = %d, s->fval = %g\n",  count, s->fval);
        
        
         if (count%200 == 0)
         {
            complx * s2 = normali2(E,0);
	  
            int test =  test_2nd_deriv(&psi_state,s2);
            freearr_(s2);

            if (test == 1)
            {
               printf("\033[32;1m Minumum found value = %lf \033[0m\n",s->fval);
               exit(2);
            }
            else if (test == -1)
               break;
         } 
         

      }
      printf(" Still not a minimum after %d iterations\n", countmax);
      printf("num_descents = %d\n", num_descents);
      printf("\033[31;1m *** \033[0m\n");
      gsl_multimin_fminimizer_free (s);
   }
   exit(2);
   gsl_vector_free(x);
   gsl_vector_free(ss);
   return 1e10;
}


complx ** makeEs(PARAMS * pm, double * evectors)
{
   int i,j,M,N;
   complx ** E;
   ull * binary_basis = enumerate_r_basis(pm->num_sites, pm->num_particles);
   ull * regions = (calc_regions(pm))[0];
   N = pm->numstates;
   M = size_(regions);
   printf("****************************************************************\n");
   printf("****************************************************************\n");
   printf("****************************************************************\n");
   printf("****************************************************************\n");
   printf("****M= %d, ***N =%d\n",M,N);
   printf("****************************************************************\n");
   new2darr_(E,N,M);
   for(i=0; i < N;i++)
   {
     double * evector = evectors+N*i;
     for(j=0; j < M;j++)
     {
        E[i][j] = evector[regions[j]];
     }
   }
   return E;
}

complx * make_coeffs(int n,int real_or_complx)
{
   int i;
   complx * ca;
   double norm;
   newarr_(ca,n);


   for(i=0;i < n;i++)
   {
      if (real_or_complx)
     {
         ca[i] = gaussian_random();
      }
      else
      {
         ca[i] = c_gaussian_random();
      }
   }
   norm = 1.0/sqrt(L2((double *) ca,2*n));

   for(i=0;i < n;i++)
      ca[i] *= norm;

   return ca;
}

double * cabs_array(complx * carray)
{
   int i;
   double * arr;
   newarr_(arr,size_(carray));
   for(i=0; i < size_(carray);i++)
      arr[i] = cabs(carray[i]);
   return arr;
}

double unitary_min(CG *cg, PARAMS * pm,  complx *coef, double * psiEs, double *
energy)
{
   complx ** E =  makeEs(pm, psiEs);
   complx ** newE;
   PSI_STATE psi_state;
   cp2darr_(newE,E);
   int N = size_(E);
   int M = size_(E[0]);
   psi_state.N = N;
   psi_state.M = M;
   psi_state.E = E;
   psi_state.newE = newE;
   psi_state.a = cabs_array(coef);
   double ** J = makeJ(&psi_state);
   psi_state.J = J;
  
 
   int size_of_box = cg->size_of_box;
   psi_state.size_of_box = size_of_box;
   int numstates = pm->numstates;
   psi_state.numstates = numstates;
   psi_state.pm = pm;
   psi_state.cg = cg;
   psi_state.psiEs = psiEs;
   psi_state.energy = energy;  
   double maxProb = newPsi(&psi_state);
   
  // double minimum = make_ran_state(N,M, E,newE, psi_state.a, J);
   return maxProb;
}


complx ** make_orth_mat(int N, int M)
{
   int i,j,k;
   complx ** E;
   new2darr_(E,N,M);
   for(i=0;i < N;i++)
     for(j=0;j < M;j++)
         E[i][j] = gaussian_random();

   
   for(j=0;j < M;j++)
   {
      for(k=0;k < j;k++)
      {
         double k_dot_j=0;
         for(i=0;i < N;i++)
            k_dot_j += E[i][j]*E[i][k];
         for(i=0;i < N;i++)
            E[i][j] -= k_dot_j*E[i][k];
      }
         double norm=0;
         for(i=0;i < N;i++)
            norm += E[i][j]*conj(E[i][j]);
         norm = 1/sqrt(norm);
         for(i=0;i < N;i++)
            E[i][j] *= norm;
   }
   //test orthnormality
   for(j=0;j < M;j++)
      for(k=0;k <= j; k++)
      {
         complx dot=0;
         for(i=0;i < N;i++)
            dot += E[i][j]*E[i][k];
         printf("E_%d dot E_%d = %lf\n",j,k,creal(dot));
      }

   return E;
}

double normE(complx ** E,int j)
{
   double norm=0;
   int i;
   for(i=0;i < size_(E);i++)
      norm += E[i][j]*conj(E[i][j]);
   return norm;
}

double maxEsite(double * a, complx ** E, int j)
{
   int i;
   double psi = 0;
   for(i=0; i < size_(E);i++)
   {
      psi += a[i]*cabs(E[i][j]);
   }
   return psi;
}

double unitary_test(int N, int M)
{

   int k;
   complx ** E =  make_orth_mat(N,M);
   complx * ca =  make_coeffs(N,1);
   double * a = cabs_array(ca);
   complx ** newE;
   freearr_(ca);
   PSI_STATE psi_state;
   cp2darr_(newE,E);
   psi_state.N = N;
   psi_state.M = M;
   psi_state.E = E;
   psi_state.newE = newE;
   psi_state.a = a;
   double ** J = makeJ(&psi_state);
   psi_state.J = J;
   double psi=maxEsite(a,E,0);
   printf("maxEsite[E,0] = %lf\n",psi*psi);
   complx * s0 = normali2(newE,0);
  // printf("energy_complex(psi_state,s0) = %lf\n",
  // energy_complex(&psi_state,s0));
#if 0
   PARAM_UNITARY param_unitary;
   psi_state.N = N;
   psi_state.M = M;
   psi_state.E = E;
   psi_state.newE = newE;
   psi_state.a = a;
   psi_state.J = J;

   double p0[3] = {0.1, 1.3, 0.6};
   param_unitary.psi_state = &psi_state;
   gsl_vector * x = gsl_vector_alloc (3);
   for(k=0;k < 3;k++)
   {
      p0[k] = 1.5*ran();
      gsl_vector_set (x, k, p0[k]);
   }
   param_unitary.i = 0;
   param_unitary.j = 3;

   minEunitary(x, &param_unitary);

   param_unitary.i = 2;
   param_unitary.j = 3;
   minEunitary(x, &param_unitary);
#endif

#if 0
   for(k=0;k< N;k++)
   {
      int j;
      s0[k] = 1;
      for(j=0;j<N;j++)
      {
         if (abs(j-k) < 2)
          J[j][k]  = 1;
         else
          J[j][k]  = 0;
      }
   }
   double * evs =   scndDerivEV(s0, J);
#endif

  newPsi(&psi_state);
  double minimum = make_ran_state(N,M, E,newE, a, J);
//  return minimum;
  return 0;
}
