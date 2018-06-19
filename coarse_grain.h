#ifndef COARSE_GRAIN_H
#define COARSE_GRAIN_H

typedef struct {
   //x coarse_graining:
   //length of system, size of each box:
   int L, size_of_box;
   int ** c_g; //c_g is the coarse graining = set of all projectors
              //cg[i][k]: i indexes the projector in x, k runs over the basis vectors that are
              //contained in that projector. The projector sig_arr[i] is the signature e.g. 2011.
              //2011 could be the state 11001010 or 11000101 etc  for size_of_box=2.
   int numCoarseEs;//Number of energy coarse graining bins
   double DeltaE;//Size of energy bucket
   int * Eindices;
   int ** sig_arr;
   ull * states;
   double ** density;//The coarse grained density
   double ** Ps;//Probablities coarse grained into bins of size_of_box,DeltaE;
} CG;


void  makeEindices(CG * cg, int numCoarseEs, int n, double * evalues);
int ** create_x_coarse_graining(CG * cg,  int num_particles, int size_of_box, ull * states);
CG * create_CG(PARAMS * pm, int size_of_box, int numCoarseEs);

#endif
