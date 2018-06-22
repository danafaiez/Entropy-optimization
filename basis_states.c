#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <complex.h>
#include "defs.h"
#include "funcs.h"
#include "data_structures.h"
#include "htable.h"
#include "coarse_grain.h"
#include "ran.h"

HTABLE * create_hash(ull * states)
{
    int numstates = size_(states);
    int i;
    
    HTABLE * ht = htable_create(numstates,sizeof(states[0]));
    for(i=0;i < numstates;i++)
    {
     ull s =  states[i];
      _htable_set(ht, &s, i);
    }
    return ht;
}

double L2(double * arr, int length)
{
  double norm = 0.0;
  int i;
  for(i=0; i < length; i++)
    norm += arr[i]*arr[i];
  return norm;
}

long long int choose (int N, int Ne)
{
   int i;
   long long int count1=1,count2=1;
   if (N < Ne)
      return 0;
   if (N == Ne)
      return 1;

   for(i=1;i<=Ne;i++)
   {
      count1 *= i;
      count2 *= N-i+1;
   }
   return(count2/count1);
}

void  print_binary_(FILE * file, int n, int numsites)
{
//   char  * c;
//   alloca_(c,64);
   int i=0;
   int * tmp;
   newarr_(tmp,numsites);
   for(i=0;i < numsites;i++) 
   {
      tmp[i] = n&1;
      n >>= 1;
   }
//   c[i] = 0;
   i=numsites;
   for(--i;i >=0; i--) 
   {
      fprintf(file, " %d",tmp[i]);
 //     c[i] = '0'+tmp[i];
   }
   freearr_(tmp);
}

void  print_binary(int n, int numsites)
{
   return print_binary_(stdout, n, numsites);
}

int iterate_r_basis(int num_sites, int num_particles, int pos[])
{
   int i,j;
   i = num_particles;

   if (pos[0]  == num_sites- num_particles)
      return 0;//done iterating

   while (pos[i-1] == pos[i] - 1)
   {
      i--;
   }
   i--;
   pos[i]++;

   for(j=i+1; j < num_particles; j++)
   {
      pos[j] = pos[i] + j-i;
   }

   return 1;
}

void high_to_low(int num_sites, int num_particles, int * pos, int * reordered)
{
   int i;
   for(i=0; i < num_particles; i++)
      reordered[i] = num_sites-1-pos[i];
}

ull positions_to_binary(int * pos, int num_particles)
{
   unsigned int binary = 0;
   int i;

   for(i=0; i < num_particles;i++) 
   {
      binary |= 1 << pos[i];
   }
   return binary;
}

ull * enumerate_r_basis(int num_sites, int num_particles)
{
   long tot=0;
   int i;
   int tot_size;
   int pos[num_particles+1];
   int reordered[num_particles];
   ull binary;
   ull * r_basis;
   int print = 0;

   //init
   for(i=0; i < num_particles; i++)
      pos[i] = i;
   pos[num_particles] = num_sites;
   tot_size = choose(num_sites,num_particles);
   if (tot_size == 0)
      return 0;

   newarr_(r_basis,tot_size);

   high_to_low(num_sites, num_particles, pos, reordered);
   binary = positions_to_binary(reordered, num_particles);

   r_basis[tot++] = binary;

   while(iterate_r_basis(num_sites,num_particles,pos))
   {
      high_to_low(num_sites, num_particles, pos, reordered);
      binary = positions_to_binary(reordered, num_particles);

      if (print)
      {
	 for(i=0; i < num_particles;i++) 
	    printf("%d ",reordered[i]);
	 print_binary(binary, num_sites);
	 printf(" binary: %lld",binary);   
	 printf(" masked: %lld", binary&((1<<5)-1));
	 printf("\n");
      }
      //assert(binary_old > binary);
      //binary_old = binary;

      r_basis[tot] = binary;
      tot++;
   }
//   printf("total = %ld\n",tot);
   return r_basis;
}

//basis1 smaller than basis2
_Complex double * convert1to2(PARAMS * pm, _Complex double * psi1,int Linit, ull * basis1, int Lfinal, ull * basis2)
{
   int i;
   int n2 = size_(basis2);
   _Complex double * psi2;
   newarr_(psi2,n2);
   for(i=0; i < size_(basis2);i++)
   {
      int numparts=  __builtin_popcountll(basis2[i]&mask(Linit));//count num of 1's in [0,Linit]
	 if (numparts == pm->num_particles)

	 {
	    int j;
	    for(j=0; j < size_(basis1);j++)
	    {
	       if (basis1[j] == basis2[i])
	       {
		  psi2[i] = psi1[j];
	       }
	    }

	 }
   }
   return psi2;
}
 

  
double scalar_product(double * psi1, double * psi2)
{
   int n = size_(psi1);
   int m = size_(psi2);
   int i=0;
   assert(n == m);
   double prod = 0.0;
   while(i < n)
   {
     prod += psi1[i]*psi2[i]; 
     i++;
   }
   return prod;
}


_Complex double scalar_product_complex_double(_Complex double * psi1, int n, double * psi2)
{
   int i=0;
   _Complex double prod = 0.0;
   while(i < n)
   {
     prod += psi1[i]*psi2[i]; 
     i++;
   }
   return prod;
}
#if 0
_Complex double * psi_thermal(int  numstates, double beta, double * energy, double * evectors, double * e)
{
   _Complex double * psi;
   double norm = 0;
   double E = 0;
   double p;
   newarr_(psi,numstates);
   int i,j;
   for(i=0;i < numstates;i++)
   {
      
     double w = exp(-0.5*beta*energy[i]); 
     double * evector = evectors+numstates*i;
      _Complex double coef;
          do {
             double rand_r = 2*ran()-1;
             double rand_i = 2*ran()-1;
             coef = rand_r + 1.0j*rand_i;   
             p = coef*conj(coef);
         } while (p >= 1.0);
      coef *= w;
      for(j=0;j < numstates;j++)
      {
         psi[j] += coef*evector[j];
      }
      E += p*energy[i];
      norm += p;
      
   }
   *e = E/norm;
   norm = 1.0/sqrt(L2((double *) psi, 2*numstates));
   for(j=0;j < numstates;j++)
   {
      psi[j] *= norm;
   }
   return psi;
}
#else
_Complex double * psi_thermal(int  numstates, double beta, double * energy, double * evectors, double * e)
{
   _Complex double * psi;
   double norm = 0;
   double E = 0;
   double p;
   newarr_(psi,numstates);
   int i,j;
   for(i=0;i < numstates;i++)
   {
      
     double w = exp(-0.5*beta*energy[i]); 
     double * evector = evectors+numstates*i;
     _Complex double coef = c_gaussian_random();
     // double coef = gaussian_random();
      coef *= w;
      for(j=0;j < numstates;j++)
      {
         psi[j] += coef*evector[j];
      }
      E += p*energy[i];
      norm += p;
      
   }
   *e = E/norm;
   norm = 1.0/sqrt(L2((double *) psi, 2*numstates));
   for(j=0;j < numstates;j++)
   {
      psi[j] *= norm;
   }
   return psi;
}
#endif


_Complex double * coeff(_Complex double * psi, int n, double * h)
{
   int i;
   _Complex double * c;
   newarr_(c,n);
   for(i=0;i<n;i++)
   {
      c[i] = scalar_product_complex_double(psi, n, h+n*i);
   }
   return c;
}

_Complex double * convert_double2complex(double * in, int N)
{
   _Complex double * out;
   int i;
   newarr_(out, N);
   for(i=0;i < N;i++)
   {
      out[i] = creal(in[i]);
   }
   return out;
}

_Complex double * psit(_Complex double * coef, double * h, double * evalues, double t)
{
   int i;
   int n = size_(coef);
   _Complex double * psi;
   newarr_(psi,n);
   for(i=0;i<n;i++)
   {
      int base = n*i;
      int j;
      _Complex double exp_iet = cexp((-1.0*I)*evalues[i]*t);
      for(j=0;j<n;j++)
      {
         psi[j] += coef[i]*h[j+base]*exp_iet;
      }
   }
   return psi;
}

int binary_to_1d_positions(int * one_d_positions, ull binary)
{
        int part_num=0;
        int pos=0;
	while(binary)
	{
	   if (binary & 1)
	   {
	      one_d_positions[part_num] = pos;
	      part_num++;
	   }
	   binary = binary >> 1;
           pos++;
	}
	return part_num;
}

void one_d_to_2d(int  Lx, COORD *c, int * oned_position, int n)
{
   int i;
   for(i=0; i < n; i++)
   {
      c[i].y = oned_position[i]/Lx;
      c[i].x =  oned_position[i] - Lx*c[i].y;
   }
}


#if 1
int num_ones(long n)//compute number of ones in binary representation of an integer, e.g. num_ones(33) = 2.
{
   return   __builtin_popcountl(n);//count num of 1's 
}
#else
int num_ones(long n)//Only use if one below isn't functional, i.e. not using gcc
{
   int tot=0;
   for(;n;n >>= 1)
      tot += n&1;
   return tot;
}
#endif


int num_ones_in_range(int start,int end, long n)//compute number of ones in binary representation of an integer, e.g. num_ones(33) = 2,
                                       // starting at begin, and ending at end.
{
   n >>= start;
   int msk = (1 << (end-start)) - 1;//mask off bits we're not interested in
   n &= msk;
   return num_ones(n);
}

ull * init_bases(int num_bath_sites, int num_sites, int num_particles)
{
   int reduced_num_sites = num_sites- num_bath_sites;
   int num_bath_particles_max;
   int reduced_num_particles_min;
   int reduced_num_particles_max;
   int num_bath_particles;
   ull * basis_size;

   num_bath_particles_max  = MIN(num_bath_sites, num_particles);
   calloc_(basis_size,num_bath_particles_max+1);

   for(num_bath_particles=0; num_bath_particles <= num_bath_particles_max; num_bath_particles++)
   {

      long num_states = choose(reduced_num_sites,num_particles-num_bath_particles);
      basis_size[num_bath_particles] = num_states;

   }

   return basis_size;
}


_Complex double ** psit_projE(CG * cg, _Complex double * coef, double * psiE, double * evalues, double t)
{
   int i;
   int n = size_(coef);
   _Complex double ** psiCG;
   double Emin = evalues[0];
   double Emax = evalues[n-1];
   int ncoarse;
//   int numCoarseEs = (Emax-Emin)/cg->DeltaE;
   int numCoarseEs = cg->numCoarseEs;
   if (cg->numCoarseEs <= 0)
      ncoarse = n;
   else 
      ncoarse = cg->numCoarseEs;
   newarr_(psiCG,ncoarse);
   if (cg->Eindices == 0)
   {
        makeEindices(cg, ncoarse, n, evalues);
   }
   for(i=0; i < ncoarse;i++)
   {
      newarr_(psiCG[i],n);
   }
   for(i=0;i<n;i++)
   {
      int base = n*i;//Starting index for ith eigenstate
      int iE, j;
      if (numCoarseEs > 0)
      {
         iE = (int)((evalues[i]-Emin)/cg->DeltaE);
         if (iE  >= numCoarseEs) iE = numCoarseEs-1;//Fixes a rounding problem for last value
      }
      else
      {
         iE = i;
      }
      _Complex double exp_iet = cexp((-1.0*I)*evalues[i]*t);
      for(j=0;j<n;j++)
      {
         psiCG[iE][j] += coef[i]*psiE[j+base]*exp_iet;
      }
   }
   return psiCG;
}

double ** calc_PsXE(PARAMS * pm, CG * cg, double * psiE, double * evalues, _Complex double * coef, double t)
{
   int i,j,ncoarse;
   int num_projs_x = size_(cg->c_g);
   int n = size_(cg->states);
   double ** Ps;
   _Complex double ** psi_projE =  psit_projE(cg,  coef, psiE,  evalues, t);
   if (cg->numCoarseEs <= 0)
      ncoarse = n;
   else 
      ncoarse = cg->numCoarseEs;

   new2darr_(Ps,ncoarse,num_projs_x);
   for(j=0;j < ncoarse;j++)
   {
      _Complex double * psi = psi_projE[j];
      for(i=0;i < num_projs_x;i++)
      {
         int k;
         for(k=0; k < size_(cg->c_g[i]);k++)
         {
            int indx = cg->c_g[i][k];
            Ps[j][i] += psi[indx]*conj(psi[indx]);
         }
      }
   }
   free2darr_(psi_projE);
   return Ps;
}

double ** calc_PsEX(PARAMS * pm, CG * cg, double * evalues, double * psiEs,_Complex double * psi)
{
   int i,ncoarse;
   int iEold = -1;
   int n = size_(cg->states);
   int numCoarseEs = cg->numCoarseEs;
   int num_projs_x = size_(cg->c_g);
   double ** Tr;
   double Emin = evalues[0];
   double Emax = evalues[n-1];
   _Complex double * ampl = 0;
   newarr_(ampl,num_projs_x);
   if (cg->numCoarseEs <= 0)
      ncoarse = n;
   else 
      ncoarse = cg->numCoarseEs;
   new2darr_(Tr,ncoarse,num_projs_x);
   if (cg->Eindices == 0)
   {
      makeEindices(cg, ncoarse, n, evalues);
   }

   for(i=0;i<n;i++)
   {
      double * psiE = psiEs  +  i*n;
      int j,iE;
      //_Complex double partialTr=0;
      if (numCoarseEs <= 0)
      {
         iE = i;
      }
      else
      {
         iE = (int)((evalues[i]-Emin)/cg->DeltaE);
         if (iE  >= numCoarseEs) iE = numCoarseEs-1;//Fixes a rounding problem for last value
      }
      for(j=0;j < num_projs_x;j++)
      {
         int k;
         for(k=0; k < size_(cg->c_g[j]);k++)
         {
            int indx = cg->c_g[j][k];
            ampl[j] += psiE[indx]*conj(psi[indx]);
         }
      }
      if (iEold < iE || i == (n-1))
      
{
         iEold = iE;
         for(j=0;j < num_projs_x;j++) Tr[iEold][j] += SQR(creal(ampl[j])) + SQR(cimag(ampl[j]));
         setarr_(ampl,0);
      }
   }
   freearr_(ampl);
   return Tr;
}

double ** calc_TrXE(PARAMS * pm, CG * cg, double * evalues, double * psiEs)
{
   int i,ncoarse;
   int n = size_(cg->states);
   int numCoarseEs = cg->numCoarseEs;
   int num_projs_x = size_(cg->c_g);
   double ** Tr;
   double Emin = evalues[0];
   double Emax = evalues[n-1];

   if (cg->numCoarseEs <= 0)
      ncoarse = n;
   else 
      ncoarse = cg->numCoarseEs;

   new2darr_(Tr,ncoarse,num_projs_x);
   if (cg->Eindices == 0)
   {
       makeEindices(cg, ncoarse, n, evalues);
   }

   for(i=0;i<n;i++)
   {
      double * psiE = psiEs  +  i*n;
      int j,iE;
      if (numCoarseEs <= 0)
      {
         iE = i;
      }
      else
      {
         iE = (int)((evalues[i]-Emin)/cg->DeltaE);
         if (iE  >= numCoarseEs) iE = numCoarseEs-1;//Fixes a rounding problem for last value
      }
      for(j=0;j < num_projs_x;j++)
      {
         int k;
         for(k=0; k < size_(cg->c_g[j]);k++)
         {
            int indx = cg->c_g[j][k];
            Tr[iE][j] += SQR(psiE[indx]);
         }
      }
   }
   return Tr;
}

void json_print_density(PARAMS * pm,CG * cg,FILE * f)
{
   int i,j,ncoarse;
   int num_projs_x = size_(cg->c_g);
   int n = size_(cg->states);
   if (cg->numCoarseEs <= 0)
      ncoarse = n;
   else 
      ncoarse = cg->numCoarseEs;
   fprintf(f,"[");
   for(j=0;j < ncoarse;j++)
   {
      fprintf(f,"[");
      for(i=0;i < size_(cg->density[j]);i++)
      {
         fprintf(f,"%lf",cg->density[j][i]);
         if (i < size_(cg->density[j]) -1) fprintf(f,",");
      }
      if (j < ncoarse) 
      {
         fprintf(f,"],");
      }
      else
      {
         fprintf(f,"]");
      }

      fprintf(f,"\n");
   }
   fprintf(f,"]");
}

double ** coarse_density(PARAMS * pm, CG * cg, double ** Ps)
{
   double ** dens=0;
   int i,j,k,E_index,ncoarse;
   int n_boxes=cg->L/cg->size_of_box;
   int n = size_(cg->states);
   if (cg->numCoarseEs <= 0)
      ncoarse = n;
   else 
      ncoarse = cg->numCoarseEs;
   new2darr_(dens,ncoarse,n_boxes);
   for(E_index=0;E_index < ncoarse;E_index++)
   {
      for(i=0; i < size_(cg->sig_arr);i++)
      {
         for(j=0; j < n_boxes;j++)
         {
            dens[E_index][j] += cg->sig_arr[i][j]*Ps[E_index][i];
         }
      }
   }
   cg->density = dens;
   return dens;
}

double ObsEntropyEX_micro(PARAMS * pm, CG * cg, double * psiEs, double * evalues,double ** evectors)
{
   assert(cg->numCoarseEs <= 0);
   int i,j,k,ncoarse;
   int num_projs_x = size_(cg->c_g);
   double S_o=0.0;
   int n = size_(cg->states);
   int num_in_ensemble = size_(evectors);
   double ** Tr, *** Ps;
   newarr_(Ps,num_in_ensemble);
   Tr = calc_TrXE(pm, cg,evalues,  psiEs);
   for(k=0;k < num_in_ensemble;k++)
   {
      _Complex double * psi = convert_double2complex(evectors[k],pm->numstates);
      Ps[k] = calc_PsEX(pm, cg, evalues, psiEs, psi);
      freearr_(psi);
   }


   for(j=0;j < n;j++)
   {
      for(i=0;i < num_projs_x;i++)
      {
         double p =0;
         double tr= Tr[j][i];
         for(k=0;k < num_in_ensemble;k++)
         {
            p +=  Ps[k][j][i];
         }
         p /= num_in_ensemble;
         if (p > 0 && tr > 0)
         {
            S_o +=  -p*log(p/tr);
         }
      }
   }
   free2darr_(Tr);
   for(k=0;k < num_in_ensemble;k++)
   {
      free2darr_(Ps[k]);
   }
   freearr_(Ps);
   return S_o;
}


double ObsEntropyEX(PARAMS * pm, CG * cg, double * psiEs, double * evalues, _Complex double * psi)
{
   int i,j,ncoarse;
   int num_projs_x = size_(cg->c_g);
   double S_o=0.0;
   int n = size_(cg->states);
   double ** Tr = calc_TrXE(pm, cg,evalues,  psiEs);
   double ** Ps = calc_PsEX(pm, cg, evalues, psiEs, psi);
   double ** density = coarse_density(pm, cg, Ps);
   cg->Ps = Ps;
   if (cg->numCoarseEs <= 0)
      ncoarse = n;
   else 
      ncoarse = cg->numCoarseEs;
   for(j=0;j < ncoarse;j++)
   {
      for(i=0;i < num_projs_x;i++)
      {
         double p =  Ps[j][i];
         if (p > 0)
         {
            S_o +=  -p*log(p/Tr[j][i]);
         }
      }
   }
   free2darr_(Tr);
//   not freed here because it is sometimes used later but it should eventually be
//   freed or cg->Ps and cg->density should be.
//   free2darr_(Ps);
//   free2darr_(density);
   return S_o;
}


double ObsEntropyXE(PARAMS * pm, CG * cg, double * psiE, double * evalues, _Complex double * coef, double t)
{
   int i,j,ncoarse;
   int num_projs_x = size_(cg->c_g);
   int n = size_(cg->states);
   double S_o=0.0;
   double ** Tr = calc_TrXE(pm, cg,evalues,  psiE);
   double ** Ps = calc_PsXE(pm, cg, psiE, evalues, coef, t);
   double ** density = coarse_density(pm, cg, Ps);
   if (cg->numCoarseEs <= 0)
      ncoarse = n;
   else 
      ncoarse = cg->numCoarseEs;
   cg->Ps = Ps;
   for(j=0;j < ncoarse;j++)
   {
      for(i=0;i < num_projs_x;i++)
      {
         double p =  Ps[j][i];
         if (p > 0)
         {
            S_o +=  -p*log(p/Tr[j][i]);
         }
      }
   }
   free2darr_(Tr);
//   free2darr_(Ps);
   return S_o;
}

ull ** calc_regions(PARAMS * pm)//for the moment, the special case of just two separate regions in Hilbert space
{
   int i;
   ull ** regions;
   unsigned long * region_sizes;
   ull * basis_size =  init_bases(pm->num_bath_sites, pm->num_sites, pm->num_particles);
   unsigned long long * binary_basis = enumerate_r_basis(pm->num_sites, pm->num_particles);
   int num_in=0, num_out=0;
   int IN=0, OUT=1;
   int x_begin= (pm->num_sites - pm->num_bath_sites)/2;
   newarr_(regions,2);
//   newarr_(regions[IN],1);
//   newarr_(regions[OUT],1);

   for(i=0; i < pm->numstates; i++)
   {
      unsigned long s = binary_basis[i];
     // if (num_ones_in_range(0, pm->num_bath_sites, s) == pm->num_particles)
      if (num_ones_in_range(x_begin, x_begin+pm->num_bath_sites, s) == pm->num_particles)
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



