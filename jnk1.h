#ifndef ENERGY_GRAINING_H
#define ENERGY_GRAINING_H

typedef struct {
   int num_parts, L, size_of_box, rank;
   double * evectors, * energy, *num;
   ull ** cellstates;
} EG;
 
EG * energy_cell_evectors(PARAMS * pm, int size_of_box);
#endif // ENERGY_GRAINING_H
