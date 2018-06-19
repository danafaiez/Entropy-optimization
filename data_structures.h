#ifndef DATA_STRUCTURES_H
#define DATA_STRUCTURES_H

#include "htable.h"

typedef struct {
   int num_sites, num_particles,
       L, Linit,x_init, x_fin, numstates, rbasis,
       kx,ky,
       size_of_box,
       num_evs,//number of eigenvectors to process when converting to real space. We do this because of RAM limitations. 
       num_bath_sites,
       num_in_ensemble,//when calculating the microcanonical density matrix, how many separate energies to include
       num_of_Ss,//the number of separate entropies that are calculated
       max_bath_sites;//how many bath sites to go to when calculating S_ent
       double t[2], U,Up;
       int minimize_S_EX,minimize_Sent,minimize_foe;
} PARAMS;

typedef struct{
   int x,y;
} COORD;

typedef struct{
   double x,y;
} COORD_DOUBLE;



#endif//DATA_STRUCTURES_H
