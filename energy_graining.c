#include <string.h>
#include <complex.h>
#include <assert.h>
#include <math.h>
#include <float.h>
#include "defs.h"
#include "io.h"
#include "fock.h"
#include "ham.h"
#include "data_structures.h"
#include "basis_states.h"
#include "funcs.h"
#include "clapack.h"
#include "ent_entropy.h"
#include "coarse_grain.h"
#include "energy_graining.h"
//#define  MAXFLOAT 1e10

double * make_product_E_state(EG * eg, int num_parts1, int ev_num1, int num_parts2, int ev_num2)
{
   int i, j;
   int size_of_box = eg->size_of_box;
   int total_box_size = 2*size_of_box;
   int tot_num_parts = num_parts1+num_parts2;
   int rank1 = size_(eg->cell_basis_states[num_parts1]);
   int rank2 = size_(eg->cell_basis_states[num_parts2]);
   double * psi1 = eg->evectors[num_parts1] +  rank1*ev_num1;
   double * psi2 = eg->evectors[num_parts2] +  rank2*ev_num2;
   ull * states = enumerate_r_basis(2*size_of_box,tot_num_parts);//only consider the case with product of two psi's
   double * psi=0;
   newarr_(psi,size_(states));
   for(i=0; i < size_(states);i++)
   {
      int numparts=  __builtin_popcountll(states[i]&mask(size_of_box));//count num of 1's in [0,total_box_size]
      if (numparts  != num_parts1)
         continue;
      int state1 = mask(size_of_box)&states[i];
      int state2 =   states[i] >> size_of_box;
      int indx1 = eg->index_n_state[num_parts1][state1];
      int indx2 = eg->index_n_state[num_parts2][state2];
      psi[i] = psi1[indx1]*psi2[indx2];
   }
   return psi;
}

EG * energy_cell_evectors(PARAMS * pm, int size_of_box)
{
      PARAMS pminit = pm[0];
      EG * eg=0;
      newarr_(eg,1);
      eg->num_parts = pm->num_particles;
      eg->L = pm->L;
      eg->size_of_box = size_of_box;
      int max_num_parts = pm->num_particles;
      int i;
      int * indices;
      double ** h=0, ** e=0;
      int * done;
      double ** evectors;
      newarr_(h,max_num_parts+1);
      newarr_(e,max_num_parts+1);
      newarr_(eg->cell_basis_states,max_num_parts+1);
      pminit.L = size_of_box;

      newarr_(eg->map2to1,max_num_parts+1);
      newarr_(eg->index_n_state,max_num_parts+1);
      for(i=0; i <= max_num_parts; i++)
      {
         int state_ind;
         eg->cell_basis_states[i] = enumerate_r_basis(size_of_box,i);
         newarr_(eg->index_n_state[i],1 << size_of_box);
         for(state_ind=0;state_ind < size_(eg->cell_basis_states[i]);state_ind++)
         {
             eg->index_n_state[i][eg->cell_basis_states[i][state_ind]] = state_ind;
         }
         int numstates = size_(eg->cell_basis_states[i]);
         eg->rank += numstates;
         newarr_(eg->map2to1[i],numstates);
         h[i] = energy(eg->cell_basis_states[i], &pminit);
         e[i] = diag(h[i],numstates,&pminit);
      }
      evectors = h;
      newarr_(indices,max_num_parts+1);
      newarr_(done,max_num_parts+1);
      int master_index=0;
      newarr_(eg->all_cell_basis_states,eg->rank);
      newarr_(eg->all_num_parts,eg->rank);
      newarr_(eg->all_evectors,eg->rank);
      newarr_(eg->all_energy,eg->rank);
      newarr_(eg->all_state_to_index, 1 << eg->L);
      while (1)
      {
	 double next_e=DBL_MAX;
	 int next_i = -1;


	 for(i=0; i <= max_num_parts; i++)
	 {
	    if (!done[i] && e[i][indices[i]] < next_e)
	    {
	       next_e = e[i][indices[i]];
	       next_i = i;
	    }
	 }
	 if (next_i < 0)
	    break;

	 eg->all_num_parts[master_index] = next_i;
	 eg->all_cell_basis_states[master_index] = eg->cell_basis_states[next_i][indices[next_i]];
	 eg->map2to1[next_i][indices[next_i]] = master_index;
	 eg->all_energy[master_index] = next_e;
	 eg->all_state_to_index[eg->all_cell_basis_states[master_index]] = master_index;

	 indices[next_i]++;
	 master_index++;
	 if (master_index >= eg->rank)
	    break;
	 if (indices[next_i] >= size_(e[next_i]))
	    done[next_i] = 1;

      }

      int n;
      for(n=0; n <= eg->num_parts;n++)
      {
         int n_states = size_(eg->cell_basis_states[n]);
         int ind;
         for(ind=0;ind < n_states;ind++)
         {
            int e_ind = eg->map2to1[n][ind];
            eg->all_evectors[e_ind] = evectors[n]+ind*n_states;
         }
      }

      eg->evectors = evectors;
      eg->energy = e;

      return eg;
}

int e_index(int dim, int length, int * xs)
{
   int i,ind=0;
   for(i=0; i < dim; i++)
   {
      ind *= length;
      ind += xs[i];
   }
   return ind;
}

int  iterate_over_es(int dim, int size, int * e_indices)
{
   int i;
   int not_done = 0;
   for(i=0; i < dim; i++)
   {
      if (e_indices[i] == (size -1))
      {
         e_indices[i] = 0;
      }
      else
      {
         e_indices[i]++;
         return 1;
      }
      not_done += e_indices[i];
   }
   return 0;//if 0, then it's done
}

int total_parts(EG * eg, int * e_indices, int num)
{
   int i;
   int tot=0;
   for(i=0;i < num; i++)
   {
      tot += eg->all_num_parts[e_indices[i]];
   }
   return tot;
}


void free_psi_Es(_Complex double ** psi_en_x)
{
   int i;
      for(i=0;i < size_(psi_en_x);i++)
            freearr_(psi_en_x[i]);
      freearr_(psi_en_x);
}

_Complex double ** transform_step_pos_to_energy(EG * eg,_Complex double * psi_en_x, ull * states, ull * left_over_state_to_index, int L)
{
   _Complex double ** psi_new_x;
   int i,j;
   assert((L % eg->size_of_box) == 0);
   int num_parts;
   int size_psi_new_x=0;
   int * left_over_size=0;//number of states in untransformed part, as a function of number of particles in that part.
   newarr_(psi_new_x,eg->rank);
   newarr_(left_over_size,eg->num_parts+1);

   for(num_parts=0;num_parts <= eg->num_parts;num_parts++)
   {
      left_over_size[num_parts] = choose(L-eg->size_of_box,num_parts);
      size_psi_new_x += left_over_size[num_parts];
   }


   for(i=0;i < eg->rank;i++)
   {
      int n_parts = eg->all_num_parts[i];//n_parts in size_of_box cell
      int n_parts_leftover = eg->num_parts - eg->all_num_parts[i];//n_parts in leftover part
      int n_states = left_over_size[n_parts_leftover];//n_states in leftover part
      newarr_(psi_new_x[i],size_psi_new_x);//psi_new_x[energy state]
   }
   freearr_(left_over_size);

   for(j=0; j < size_(states);j++)
   {
      ull state = states[j];
      int cell_state = state >> (L-eg->size_of_box);
      int total_cell_index = eg->all_state_to_index[cell_state];
      int left_over_state = state & mask(L-eg->size_of_box);
      int left_over_index = left_over_state_to_index[left_over_state];
      int num_parts_cell = eg->all_num_parts[total_cell_index];
      int n_states = size_(eg->cell_basis_states[num_parts_cell]);
      int all_vect_ind =  eg->index_n_state[num_parts_cell][cell_state];
      //   Compute psi_new_x for all possible energies. Impossible e_ind's would have the wrong number
      //   of particles in the system. What you'd want is a for loop over all e_ind:
      //   psi_new_x[e_ind][left_over_index] += psi_en_x[j]*eg->evectors_aranged_in_order[e_ind][total_cell_index];
      //   but that'd be inefficient because evectors_aranged_in_order would be zero for e_inds that
      //   have the wrong number of particles, i.e. numparts[total_cell_index] + numparts[e_ind] !=    numparts_total
      for(i=0; i < n_states ;i++)
      {
         int e_ind = eg->map2to1[num_parts_cell][i]; 
         psi_new_x[e_ind][left_over_index] += psi_en_x[j]*eg->all_evectors[e_ind][all_vect_ind];
      }

   }
   return psi_new_x;
}

_Complex double * transform_pos_to_energy(EG * eg, _Complex double * psit)
{
   int num_x_bins = eg->L/eg->size_of_box;
   int * e_indices;
   int i;
   _Complex double ** psi_en_x, ** psi_en_x_next, ** psi_en_x_new;
   newarr_(e_indices,num_x_bins);
   ull * states_left_over_n;
   ull * states_left_over=0;
   ull * states = enumerate_r_basis(eg->L,eg->num_parts);
   newarr_(psi_en_x,1);
   memcpyarr_(psi_en_x[0],psit);
   int x_bin;

   for(x_bin=0;x_bin < num_x_bins; x_bin++)
   {
      ull * left_over_state_to_index;
      int shift = x_bin*eg->size_of_box;
      newarr_(left_over_state_to_index,1 << (eg->L-shift-eg->size_of_box));
      int num;
      int count=0;
      setarr_(e_indices,0);
      newarr_(psi_en_x_new,pow(eg->rank,x_bin+1));
      for(num=0;num <= eg->num_parts;num++)
      {
         int state_index;
         states_left_over_n = enumerate_r_basis(eg->L-shift-eg->size_of_box,num);
         if (states_left_over_n)
         {
            for(state_index=0;state_index < size_(states_left_over_n);state_index++)
            {
               ull binary_state = states_left_over_n[state_index];
               left_over_state_to_index[binary_state] = count;
               appendarr_(states_left_over, binary_state);
               count++;
            }
         }
         freearr_(states_left_over_n);
      }

      do
      {
         int ind = e_index(x_bin,eg->rank,e_indices);
         //num = eg->num_parts-total_parts(eg, e_indices,x_bin);
         //psi_en_x_next[ind] is the evector
         //transform the x-space piece of the evector psi_en_x[ind] at a given set of energy indices
         //(ind). Input the reduces basis states states_left_over, and num of parts left over, and
         //length. Output will have transformed the next size_of_box sites to the local energy
         //basis.
         psi_en_x_next = transform_step_pos_to_energy(eg,psi_en_x[ind],states,left_over_state_to_index, eg->L-shift);
         for(ind=0;ind < eg->rank;ind++)
         {
            e_indices[x_bin] = ind;
            ull combined_index = e_index(x_bin+1,eg->rank,e_indices);
            psi_en_x_new[combined_index] = psi_en_x_next[ind];//larger e part and smaller x part then psi_en_x.
         }
	 freearr_(psi_en_x_next);

      } while (iterate_over_es(x_bin, eg->rank, e_indices));

      free_psi_Es(psi_en_x);
      psi_en_x = psi_en_x_new;

      freearr_(states);
      states = states_left_over;
      states_left_over = 0;
      freearr_(left_over_state_to_index);
   }
   freearr_(states);
   freearr_(e_indices);
   _Complex double * psi=0;
   newarr_(psi,size_(psi_en_x_new));
   for(i=0;i < size_(psi_en_x_new);i++)
      psi[i] = psi_en_x_new[i][0];
   free2darr_(psi_en_x_new);
   return psi;
}

double Sobs_fine_grain_E(_Complex double * psiE)
{
   int i;
   double S=0;
   for(i=0; i < size_(psiE);i++)
   {
      double absPsi_r = creal(psiE[i]);
      double absPsi_i = cimag(psiE[i]);
      double p = SQR(absPsi_r) + SQR(absPsi_i);
      if (p > 0)
         S += p*log(p);
   }
   return -S;
}

double Sobs_fine_grain_micro_E(_Complex double ** psiE)
{
   int i;
   double S=0;
   int num_evs = size_(psiE);
   int length = size_(psiE[0]);
   double norm = 1.0/num_evs;
   for(i=0; i < length;i++)
   {
      double p = 0;
      int j;
      for(j=0;j < num_evs;j++)
      {
         double absPsi_r = creal(psiE[j][i]);
         double absPsi_i = cimag(psiE[j][i]);
         p += SQR(absPsi_r) + SQR(absPsi_i);
      }
      p *= norm;
      if (p > 0)
         S += p*log(p);
   }
   return -S;
}

_Complex double * superpose_psi(double ** psiE)
{
   int i,j,k;
   double S;
   int num_evs = size_(psiE);
   int length = size_(psiE[0]);
   double norm = 0;
   double p;
   _Complex double * psi;
   newarr_(psi,length);
   for(j=0;j < num_evs;j++)
   {
      _Complex double coef;
      do 
      {
         double rand_r = 2*ran()-1;
         double rand_i = 2*ran()-1;
         coef = rand_r + 1.0j*rand_i;
         p = coef*conj(coef);
      } while (p >= 1.0);

      for(k=0;k < length;k++)
      {
         psi[k] += coef*psiE[j][k];
      }
//      norm += p;
   }

   norm = 1.0/sqrt(L2((double *) psi, 2*length));
   //norm = 1.0/sqrt(norm);
   for(j=0;j < length;j++)
   {
      psi[j] *= norm;
   }

   return psi;
}

