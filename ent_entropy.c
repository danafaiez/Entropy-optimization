#if 1

#include <math.h>
#include <assert.h>
#include <complex.h>
#include "basis_states.h"
#include "data_structures.h"
#include "defs.h"
//#include "io.h"

long indx(long y_size,long i, long j)
{
   return j*y_size + i;//fortran convention
}

_Complex double ** callocate_rho_cmplx(int num_bath_particles_max, ull * basis_size)
{
   _Complex double ** rho;
   int i;
   calloc_(rho, num_bath_particles_max+1);
   for(i=0; i < num_bath_particles_max+1; i++)
      calloc_(rho[i], ((basis_size[i]+1)*(basis_size[i]))/2);//Since it's Hermitian we only store 1/2 the matrix.

   return rho;
}


double ** callocate_rho(int num_bath_particles_max, ull * basis_size)
{
   double ** rho;
   int i;
   calloc_(rho, num_bath_particles_max+1);
   for(i=0; i < num_bath_particles_max+1; i++)
      calloc_(rho[i], ((basis_size[i]+1)*(basis_size[i]))/2);//Since it's Hermitian we only store 1/2 the matrix.

   return rho;
}

_Complex double ** red_dens_mat_complex(int num_sites, int num_particles, int num_bath_sites, ull * binary_basis, unsigned  long nbasis, _Complex double * evector) 
{

   _Complex double ** rho;
   int num_bath_particles_max  = MIN(num_bath_sites, num_particles);
   int i;
   unsigned long mask;
   unsigned long base=0;
   ull * basis_size =  init_bases(num_bath_sites, num_sites, num_particles);

   rho = callocate_rho_cmplx(num_bath_particles_max, basis_size);
#if 0
   calloc_(rho, num_bath_particles_max+1);
   for(i=0; i < num_bath_particles_max+1; i++)
      calloc_(rho[i], ((basis_size[i]+1)*(basis_size[i]))/2);//Since it's Hermitian we only store 1/2 the matrix.
#endif

   //create a mask that's only looks at the number of particles in the bath region
   // i.e. 11111100000000
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
	    //rho[nbath][indx(basis_size[nbath],i,j)] += evector[base+i]*conj(evector[base+j]);
	    rho[nbath][count] += evector[base+i]*conj(evector[base+j]);
	    //rho[nbath][count] += ran();
	    count++;
	 }
      }
      base += basis_size[nbath];
   }
   free(basis_size);
   return rho;
}


double ** red_dens_mat(int num_sites, int num_particles, int num_bath_sites, ull * binary_basis, unsigned  long numstates, double * evector) 
{

   double ** rho;
   int num_bath_particles_max  = MIN(num_bath_sites, num_particles);
   int i;
   unsigned long mask;
   unsigned long base=0;
   ull * basis_size =  init_bases(num_bath_sites, num_sites, num_particles);

   rho = callocate_rho(num_bath_particles_max, basis_size);
#if 0
   calloc_(rho, num_bath_particles_max+1);
   for(i=0; i < num_bath_particles_max+1; i++)
      calloc_(rho[i], ((basis_size[i]+1)*(basis_size[i]))/2);//Since it's Hermitian we only store 1/2 the matrix.
#endif

   //create a mask that's only looks at the number of particles in the bath region
   // i.e. 11111100000000
   mask  = ((1 << num_bath_sites) -1) << (num_sites-num_bath_sites);


   while (base < numstates)
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
	    //rho[nbath][indx(basis_size[nbath],i,j)] += evector[base+i]*conj(evector[base+j]);
	    rho[nbath][count] += evector[base+i]*evector[base+j];
	    //rho[nbath][count] += ran();
	    count++;
	 }
      }
      base += basis_size[nbath];
   }
   free(basis_size);
   return rho;
}

double * read_evectors(char * evector_in_name, PARAMS * pm);

int diag_rho(int N, double * matrix_double, double * eigenvalues);
int diag_rho_cmplx(int n, double * matrix_double, double * eigenvalues);

double * calc_ent_entropy_(double * evectors, PARAMS * pm, int num_bath_sites, int num_in_ensemble)
{
   ull * basis_size =  init_bases(pm->num_bath_sites, pm->num_sites, pm->num_particles);
   ull * binary_basis = enumerate_r_basis(pm->num_sites, pm->num_particles);
   int num_bath_particles;
   int num_bath_particles_max;
   int ev_index;//index of the eigenvector
   int num_of_Ss = pm->numstates/num_in_ensemble;
   double * S;

   pm->num_bath_sites = num_bath_sites;// May want to specify a different number // bath sites

   calloc_(S,pm->numstates);

   num_bath_particles_max  = MIN(pm->num_bath_sites, pm->num_particles);


   int S_index;
   for(S_index=0; S_index < num_of_Ss; S_index++)
   {
      int j;
      double ** rho_micro = callocate_rho(num_bath_particles_max, basis_size);

      for(j=0;  j < num_in_ensemble; j++)
      {
	 int i;
	 int ev_index = S_index*num_in_ensemble + j;

	 //to start from random superposition of neighboring energy
	 //eigenvectors, put a function call here to add evectors together and
	 //replace the last argument of red_dens_mat:
	 double ** red_rho = red_dens_mat(pm->num_sites, pm->num_particles, pm->num_bath_sites, binary_basis, pm->numstates,  evectors + ev_index*pm->numstates);

	 for(num_bath_particles=0; num_bath_particles < num_bath_particles_max+1; num_bath_particles++)
	 {
	    int count;
	    for(count = 0; count < basis_size[num_bath_particles]*(basis_size[num_bath_particles]+1)/2; count++)
	    {
	       rho_micro[num_bath_particles][count] += red_rho[num_bath_particles][count];
	    }
	 }
	 for(i=0; i < num_bath_particles_max+1; i++)
	    free(red_rho[i]);
      }
      for(num_bath_particles=0; num_bath_particles < num_bath_particles_max+1; num_bath_particles++)
      {
	 int count;
	 for(count = 0; count < basis_size[num_bath_particles]*(basis_size[num_bath_particles]+1)/2; count++)
	    rho_micro[num_bath_particles][count] /= num_in_ensemble;
      }

      for(num_bath_particles=0; num_bath_particles < num_bath_particles_max+1; num_bath_particles++)
      {
	 double * eigenvalues;
	 double S_=0;
	 int i, n = basis_size[num_bath_particles];
	 calloc_(eigenvalues, n);
	 diag_rho(n, (double *) rho_micro[num_bath_particles], eigenvalues);

	 for(i=0; i < n; i++)
	 {
	    assert(eigenvalues[i] > -1.0e-10);
	    if (eigenvalues[i] > 0)
	       S_ += -eigenvalues[i] * log(eigenvalues[i]);
	 }

	 S[S_index] += S_;
	 free(eigenvalues);
      }
      int i;
      for(i=0; i < num_bath_particles_max+1; i++)
	 free(rho_micro[i]);
   }
   return S;
}

double calc_ent_entropy_one_ev_complex_(_Complex double * evector, PARAMS * pm, int num_bath_sites)
{
   ull * basis_size =  init_bases(pm->num_bath_sites, pm->num_sites, pm->num_particles);
   ull * binary_basis = enumerate_r_basis(pm->num_sites, pm->num_particles);
   int num_bath_particles;
   int num_bath_particles_max;
   int ev_index;//index of the eigenvector
   double  S=0;
   int i;

   pm->num_bath_sites = num_bath_sites;// May want to specify a different number // bath sites

   num_bath_particles_max  = MIN(pm->num_bath_sites, pm->num_particles);


   {
      _Complex double ** red_rho = red_dens_mat_complex(pm->num_sites, pm->num_particles, pm->num_bath_sites, binary_basis, pm->numstates,  evector);
      
      for(num_bath_particles=0; num_bath_particles < num_bath_particles_max+1; num_bath_particles++)
      {
	 double * eigenvalues;
	 double S_=0;
	 int i, n = basis_size[num_bath_particles];
	 //printf("n=basis_size[num_bath_particles]=%d\n",n);
	 calloc_(eigenvalues, n);
	 diag_rho_cmplx(n, (double *) red_rho[num_bath_particles], eigenvalues);

	 for(i=0; i < n; i++)
	 {
	    //printf("num_bath_particles:%d;",num_bath_particles);
	    //printf("i:%d;\n",i);
	   
	    assert(eigenvalues[i] > -1.0e-10);
	    if (eigenvalues[i] > 0)
	       S_ += -eigenvalues[i] * log(eigenvalues[i]);
	       //printf("%lf,",eigenvalues[i]);
	 }
	 S += S_;
	 free(eigenvalues);
      }
      //printf("\n");
      int c;
      
      for(i=0; i < num_bath_particles_max+1; i++)
	  free(red_rho[i]);
      free(red_rho);  
   }
   //printf(",\n");
   freearr_(binary_basis);
   free(basis_size);
   
   return S;
}



double calc_ent_entropy_one_ev_(double * evector, PARAMS * pm, int num_bath_sites)
{
   ull * basis_size =  init_bases(pm->num_bath_sites, pm->num_sites, pm->num_particles);
   ull * binary_basis = enumerate_r_basis(pm->num_sites, pm->num_particles);
   int num_bath_particles;
   int num_bath_particles_max;
   int ev_index;//index of the eigenvector
   double  S=0;
   int i;

   pm->num_bath_sites = num_bath_sites;// May want to specify a different number // bath sites

   num_bath_particles_max  = MIN(pm->num_bath_sites, pm->num_particles);


   {
      double ** red_rho = red_dens_mat(pm->num_sites, pm->num_particles, pm->num_bath_sites, binary_basis, pm->numstates,  evector);
      for(num_bath_particles=0; num_bath_particles < num_bath_particles_max+1; num_bath_particles++)
      {
	 double * eigenvalues;
	 double S_=0;
	 int i, n = basis_size[num_bath_particles];
	 calloc_(eigenvalues, n);
	 diag_rho(n, (double *) red_rho[num_bath_particles], eigenvalues);

	 for(i=0; i < n; i++)
	 {
	    assert(eigenvalues[i] > -1.0e-10);
	    if (eigenvalues[i] > 0)
	       S_ += -eigenvalues[i] * log(eigenvalues[i]);
	 }

	 S += S_;
	 free(eigenvalues);
      }
      for(i=0; i < num_bath_particles_max+1; i++)
	 free(red_rho[i]);
   }

   freearr_(binary_basis);
   free(basis_size);
   return S;
}

double * calc_ent_entropy_one_ev(double * evector, PARAMS * pm)
{
   double * entropies;
   double norm=0.0;
   int i;
//   int max_bath_sites = pm->num_sites/2;
   int max_bath_sites = pm->max_bath_sites;
   int num_bath_sites;

   norm = sqrt(1.0/L2((double *)evector,  pm->numstates));
   for(i=0; i < pm->numstates; i++)
      evector[i] *= norm;

   calloc_(entropies,max_bath_sites+1);
   for(num_bath_sites = pm->num_sites; num_bath_sites >= pm->num_sites - max_bath_sites; num_bath_sites--)
   {
      pm->num_bath_sites = num_bath_sites;
      entropies[pm->num_sites - num_bath_sites] = calc_ent_entropy_one_ev_(evector, pm, num_bath_sites);
   }

   return entropies;
}

double calc_ent_entropy_n_ev_(double ** evectors, PARAMS * pm, int num_bath_sites, int n)
{
   ull * basis_size =  init_bases(pm->num_bath_sites, pm->num_sites, pm->num_particles);
   ull * binary_basis = enumerate_r_basis(pm->num_sites, pm->num_particles);
   int num_bath_particles;
   int num_bath_particles_max;
   int ev_index;//index of the eigenvector
   double  S=0;
   int i,j;

   pm->num_bath_sites = num_bath_sites;// May want to specify a different number // bath sites

   num_bath_particles_max  = MIN(pm->num_bath_sites, pm->num_particles);


   {
      double ** rho_micro = callocate_rho(num_bath_particles_max, basis_size);

      for(j=0;  j < n; j++)
      {
	 int i;
	 int ev_index = j;

         double ** red_rho = red_dens_mat(pm->num_sites, pm->num_particles, pm->num_bath_sites, binary_basis, pm->numstates,  evectors[j]);

	 for(num_bath_particles=0; num_bath_particles < num_bath_particles_max+1; num_bath_particles++)
	 {
	    int count;
	    for(count = 0; count < basis_size[num_bath_particles]*(basis_size[num_bath_particles]+1)/2; count++)
	    {
	       rho_micro[num_bath_particles][count] += red_rho[num_bath_particles][count];
	    }
	 }
	 for(i=0; i < num_bath_particles_max+1; i++)
	    free(red_rho[i]);
      }

      for(num_bath_particles=0; num_bath_particles < num_bath_particles_max+1; num_bath_particles++)
      {
	 int count;
	 for(count = 0; count < basis_size[num_bath_particles]*(basis_size[num_bath_particles]+1)/2; count++)
	    rho_micro[num_bath_particles][count] /= n;
      }

      for(num_bath_particles=0; num_bath_particles < num_bath_particles_max+1; num_bath_particles++)
      {
	 double * eigenvalues;
	 double S_=0;
	 int i, n = basis_size[num_bath_particles];
	 calloc_(eigenvalues, n);
	 diag_rho(n, (double *) rho_micro[num_bath_particles], eigenvalues);

	 for(i=0; i < n; i++)
	 {
	    assert(eigenvalues[i] > -1.0e-10);
	    if (eigenvalues[i] > 0)
	       S_ += -eigenvalues[i] * log(eigenvalues[i]);
	 }

	 S += S_;
	 free(eigenvalues);
      }
      for(i=0; i < num_bath_particles_max+1; i++)
	 free(rho_micro[i]);
   }

   freearr_(binary_basis);
   free(basis_size);
   return S;
}

double * calc_ent_entropy_n_evs(double ** evectors, PARAMS * pm, int n)
{
   double * entropies;
   double norm=0.0;
   int i;
//   int max_bath_sites = pm->num_sites/2;
   int max_bath_sites = pm->max_bath_sites;
   int num_bath_sites;
   int j;

   for(j=0; j < n; j++)
   {

      norm = sqrt(1.0/L2((double *)evectors[j],  pm->numstates));
      for(i=0; i < pm->numstates; i++)
	 evectors[j][i] *= norm;
   }

   calloc_(entropies,max_bath_sites+1);
   for(num_bath_sites = pm->num_sites; num_bath_sites >= pm->num_sites - max_bath_sites; num_bath_sites--)
   {
      pm->num_bath_sites = num_bath_sites;
      entropies[pm->num_sites - num_bath_sites] = calc_ent_entropy_n_ev_(evectors, pm, num_bath_sites,n);
   }

   return entropies;
}

double * calc_ent_entropy_n_evs_superpose(double ** evectors, PARAMS * pm, int n)
{
   double * entropies;
   double norm=0.0;
   int i,j,k;
   //   int max_bath_sites = pm->num_sites/2;
   int max_bath_sites = pm->max_bath_sites;
   int num_bath_sites;
   double * superposed_evector;
   calloc_(superposed_evector, pm->numstates);

   for(j=0; j < n; j++)
   {

      norm = sqrt(1.0/L2((double *)evectors[j],  pm->numstates));
      for(i=0; i < pm->numstates; i++)
	 evectors[j][i] *= norm;
   }

   for(j=0;  j < n; j++)
   {
      int k;
      double crand = 2*ran()-1;
      for(k=0; k < pm->numstates;k++)
	 superposed_evector[k] += crand*evectors[j][k];
   }

   norm = sqrt(1.0/L2((double *)superposed_evector, pm->numstates));
   for(k=0; k < pm->numstates; k++)
      superposed_evector[k] *= norm;

   calloc_(entropies,max_bath_sites+1);
   for(num_bath_sites = pm->num_sites; num_bath_sites >= pm->num_sites - max_bath_sites; num_bath_sites--)
   {
      pm->num_bath_sites = num_bath_sites;
      entropies[pm->num_sites - num_bath_sites] = calc_ent_entropy_n_ev_(&superposed_evector, pm, num_bath_sites,1);
   }

   free(superposed_evector);
   return entropies;
}



double ** calc_ent_entropy(double * evectors, PARAMS * pm)
{
   double ** entropies;
   double norm=0.0;
   int i;
   int max_bath_sites = pm->num_sites;
   int num_bath_sites;

   norm = sqrt(1.0/L2((double *)evectors,  pm->numstates));
   for(i=0; i < pm->numstates*pm->rbasis; i++)
       evectors[i] *= norm;

   calloc_(entropies,max_bath_sites+1);
   for(num_bath_sites = pm->num_sites; num_bath_sites >= pm->num_sites - max_bath_sites; num_bath_sites--)
     {
       pm->num_bath_sites = num_bath_sites;
      entropies[pm->num_sites - num_bath_sites] = calc_ent_entropy_(evectors, pm, num_bath_sites, pm->num_in_ensemble);
     }

   free(evectors);
   return entropies;
}

double calc_dom_entropy_region(double * state_vector, PARAMS * pm, ull * region)
{
    int i; 
    double prob = 0.0;
    for(i=0; i < size_(region); i++) 
    {
       prob += SQR(state_vector[2*region[i]])+SQR(state_vector[2*region[i]+1]);
    }
    if (prob == 0)
      return 0;
    return -prob*log(prob/size_(region));
}

double calc_dom_entropy_regions(double * state_vector, PARAMS * pm, ull ** regions)
{
   int i;
   double S=0;
   for(i=0; i < size_(regions); i++)
   {
      S += calc_dom_entropy_region(state_vector, pm, regions[i]);
   }
   return S;
}

double calc_dom_entropy_block(double * state_vector, PARAMS * pm)
{
   ull  ** regions = calc_regions(pm);
   return calc_dom_entropy_regions(state_vector, pm, regions);
}


#endif
