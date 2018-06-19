#ifndef FOCK_H
#define FOCK_H

typedef unsigned long long ull;

#define mask(i) ((1 << (i)) -1)
#define bit_val(num,pos) (1&((num) >> (pos)))
#define flip_bit(num,pos) ((num) ^= (1 << (pos)))

#endif// FOCK_H
