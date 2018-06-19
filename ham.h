#ifndef HAM_H
#define HAM_H

#define LEFT  0
#define RIGHT 1
#define H(i,j) h[(i) + numstates*(j)]

double * energy(ull * states, PARAMS * pm);

#endif// HAM_H
