#ifndef ENERGY_GRAINING_H
#define ENERGY_GRAINING_H

typedef struct {
   int num_parts, L, size_of_box, rank;
   //     Inside a cell:
   double ** evectors, ** energy;//energy[num_parts][energy_index], evectors[num_parts][energy_index][position_basis_indx]
   double ** all_evectors, * all_energy;//flatten above combining num_parts and energy_index into a single index
   int * all_num_parts;//num of particles in each state above.
   ull ** cell_basis_states;//cell_basis_states[part_num][basis states_index]
   ull *  all_cell_basis_states;//above array flattened.
   int * all_state_to_index;//binary state -> master index
   int ** map2to1;//map2to1[part_num][energy_index] -> single index
   int ** index_n_state;//index_n_state[num_parts][bin_state] -> state basis index with num_parts particles
} EG;
 
EG * energy_cell_evectors(PARAMS * pm, int size_of_box);
_Complex double * transform_pos_to_energy(EG * eg, _Complex double * psit);
void free_psi_Es(_Complex double ** psi_en_x);
double Sobs_fine_grain_E(_Complex double * psiE);
double * make_product_E_state(EG * eg, int num_parts1, int ev_num1, int num_parts2, int ev_num2);
double Sobs_fine_grain_micro_E(_Complex double ** psiE);
_Complex double * superpose_psi(double ** psiE);
#endif // ENERGY_GRAINING_H
