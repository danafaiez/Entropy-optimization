#include <stdlib.h>
#include "locate_multivs.h"

typedef struct {
   int ind;
   double val;
} IND_VAL;

int psi2_comp(const void * a, const void * b)
{
   double ad = ((IND_VAL *) a)->val;
   double bd = ((IND_VAL *) b)->val;
   double  diff =  -(ad-bd); 
   if (diff > 0) return 1; 
      else if (diff == 0) return 0;
         else return -1;
}

int * top_psi_indices(PARAMS *pm, double * evector, int top_n, int is_cmplx) 
{
    unsigned long long * binary_basis = enumerate_r_basis(pm->num_sites, pm->num_particles);
    int n = pm->numstates;
    int i;
    IND_VAL * iv;
    newarr_(iv,n);
    for(i=0; i < n; i++)
    {
       iv[i].val = is_cmplx? SQR(evector[2*i])+SQR(evector[2*i+1]) : SQR(evector[i]);
       iv[i].ind = i;
    }
    qsort(iv,n,sizeof(iv[0]), psi2_comp);
    int * result;
    newarr_(result,top_n);
    for(i=0; i < top_n; i++)
    {
        result[i] = iv[i].ind;
    }

    freearr_(iv);
    return result;
}


void hist_evector(PARAMS * pm, double * evector, int is_cmplx) 
{

    unsigned long long * binary_basis = enumerate_r_basis(pm->num_sites, pm->num_particles);
    int n = pm->numstates;
    int i;
    IND_VAL * iv;
    HIST * h =  new_hist(0.0, 15*1.0/n, 32);
    calloc_(iv,n);
    for(i=0; i < n; i++)
    {
       float rho = is_cmplx ? SQR(evector[2*i])+SQR(evector[2*i+1]) : SQR(evector[i]);
       add_hist(h,rho);
    }
    printf("numstates = %d, largest prob is %d\n",n, h->max);
    FILE * f = fopen("histout","w");
    print_hist(h, f);
    fclose(f);
}


