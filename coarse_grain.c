#include <stdio.h>
#include <math.h>
#include "defs.h"
#include "htable.h"
#include "funcs.h"
#include "fock.h"
#include "basis_states.h"
#include "coarse_grain.h"


//for each state it outputs a string that corresponds to the name of a projector
int * make_signature(int L,int size_of_box,ull ullstate)
{
   int * signature;
   int n_boxes =L/size_of_box;
   int i,j;
   newarr_(signature,n_boxes);
   for(i=0;i<n_boxes;i++)
   {
      int counter=0;
      for(j=0;j<size_of_box;j++)
      {
         counter += ullstate & 1;
         ullstate = ullstate >> 1;
      }
      signature[i]=counter;
   }
   return signature;
}

int ** create_x_coarse_graining(CG * cg,  int num_particles, int size_of_box, ull * states)
{
   int numstates = size_(states);
   int L = cg->L;
   int n_boxes=L/size_of_box;
   int n_projectors=choose(num_particles + n_boxes-1,num_particles);
   int cg_index;
   int ** c_g=0;//coarse graining consisting of projectors. 
              //cg[i][k]: i indexes the projector in x, k runs over the basis vectors that are
              //contained in that projector. The projector sig_arr[i] is the signature e.g. 2011.
              //2011 could be the state 11001010 or 11000101 etc  for size_of_box=2.
   int * projector;//projectors
   int ** sig_arr=0;
   int i;
   int j=0;
   HTABLE * ht = htable_create(n_projectors,(n_boxes*sizeof(int)));
   //create space to store the signatures of n_projectors
   //newarr_(c_g,n_projectors);
   for(i=0;i<numstates;i++)//go through all states
   {
      //for each state, create a signature to find to which projector the state belongs to
      int * signature=make_signature(L,size_of_box,states[i]);
      int * p_value = _htable_get(ht,signature);//find index of projector with given signature
      if (p_value==0)
      {
         _htable_set (ht,signature,j);
         cg_index = j;
         j++;
         appendarr_(sig_arr,signature);
         appendarr_(c_g,0);
      }
      else
      {
         cg_index=*p_value;//find index of projector with given signature
      }
      //insert index of the state (i) into the projector with an index (cg_index)
      appendarr_(c_g[cg_index],i);
}
   cg->sig_arr = sig_arr;
   htable_free(ht);
   //free2darr_(sig_arr); 
   return c_g;

}

void  makeEindices(CG * cg, int numCoarseEs, int n, double * evalues)
{
   int * ind;
   int i;
   double Emin = evalues[0];
   double Emax = evalues[n-1];
   cg->numCoarseEs = numCoarseEs;
   if (numCoarseEs <= 0) //meaning no coarse graining, just fine graining
   {
      cg->DeltaE = 0;
      cg->DeltaE = (Emax-Emin)/n;
      if (cg->Eindices == 0)
      {
         newarr_(ind,n);
         for(i=0; i < n;i++)
         {
            ind[i] = i;
         }
         cg->Eindices = ind;
      }
   }
   else
   {
      cg->DeltaE = (Emax-Emin)/numCoarseEs;
      if (cg->Eindices == 0)
      {
         newarr_(ind,numCoarseEs);
         for(i=0; i < n;i++)
         {
            int j = (evalues[i]-Emin)/cg->DeltaE;
            if (j == numCoarseEs) j--;//fixes rounding problem for last value
            ind[j]++;
         }
         cg->Eindices = ind;
      }
   }
   
}

CG * create_CG(PARAMS * pm, int size_of_box, int numCoarseEs)
{
   ull * states = enumerate_r_basis(pm->L, pm->num_particles);
   CG * cg;
   newarr_(cg,1);
   cg->L = pm->L;
   cg->size_of_box = size_of_box;
   cg->numCoarseEs = numCoarseEs;
   cg->states = states;
   cg->c_g = create_x_coarse_graining(cg, pm->num_particles, size_of_box, states);
   return cg; 
}

